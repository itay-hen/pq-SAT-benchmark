/*
 * sat_benchmark.cpp — Benchmark SAT solvers on DIMACS CNF instances.
 *
 * Runs Kissat, CaDiCaL, and Glucose on every .cnf file found in the
 * same directory as this executable. Reports solving time, result
 * (SAT/UNSAT/TIMEOUT/ERROR), and optionally verifies the solution
 * against a planted-solution file if one exists.
 *
 * Output is a tab-separated file (TSV) with one row per instance.
 * Columns: instance, num_vars, num_clauses, then for each solver
 * two columns: <solver>_time and <solver>_result.
 *
 * Key robustness fixes in this version:
 *   - Proper shell-quoting of solver paths, CNF paths, and temp files.
 *   - Safer handling of std::system return codes.
 *   - Automatic use of `timeout` or `gtimeout` if available.
 *   - Fallback to running without timeout when neither exists.
 *
 * Compile:
 *   g++ -std=c++17 -O2 -o sat_benchmark sat_benchmark.cpp
 *
 * Usage:
 *   ./sat_benchmark                          # run all solvers, 60s timeout
 *   ./sat_benchmark --timeout 300            # 5-minute timeout
 *   ./sat_benchmark --solver kissat          # run only kissat
 *   ./sat_benchmark --file instance.cnf      # run on a single file
 *   ./sat_benchmark --tsv results.tsv        # write TSV output
 */

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <sys/wait.h>
#include <unistd.h>
#include <vector>

namespace fs = std::filesystem;

// ---------------------------------------------------------------------------
// Result of a single solver run
// ---------------------------------------------------------------------------

enum class Result { SAT, UNSAT, TIMEOUT, ERROR, UNKNOWN };

const char* result_str(Result r) {
    switch (r) {
        case Result::SAT:     return "SAT";
        case Result::UNSAT:   return "UNSAT";
        case Result::TIMEOUT: return "TIMEOUT";
        case Result::ERROR:   return "ERROR";
        default:              return "UNKNOWN";
    }
}

struct SolverRun {
    std::string solver_name;
    std::string instance;
    Result result = Result::UNKNOWN;
    double elapsed_sec = 0.0;
    bool solution_verified = false;
    bool verification_attempted = false;
    bool factorization_verified = false;
    bool factorization_attempted = false;
    long long recovered_p = 0, recovered_q = 0;
};

// ---------------------------------------------------------------------------
// Small utilities
// ---------------------------------------------------------------------------

std::string shell_quote(const std::string& s) {
    std::string out = "'";
    for (char c : s) {
        if (c == '\'') {
            out += "'\\''";
        } else {
            out += c;
        }
    }
    out += "'";
    return out;
}

std::string read_text_file(const std::string& path) {
    std::ifstream in(path);
    if (!in.is_open()) return {};
    return std::string((std::istreambuf_iterator<char>(in)),
                       std::istreambuf_iterator<char>());
}

std::string trim_copy(std::string s) {
    while (!s.empty() && (s.back() == '\n' || s.back() == '\r' || s.back() == ' ' || s.back() == '\t')) {
        s.pop_back();
    }
    size_t i = 0;
    while (i < s.size() && (s[i] == ' ' || s[i] == '\t')) {
        ++i;
    }
    return s.substr(i);
}

std::string find_command_in_path(const std::string& command) {
    std::string cmd = "command -v " + shell_quote(command) + " 2>/dev/null";
    FILE* pipe = popen(cmd.c_str(), "r");
    if (!pipe) return {};

    char buf[512];
    std::string result;
    while (fgets(buf, sizeof(buf), pipe)) result += buf;
    pclose(pipe);

    result = trim_copy(result);
    if (result.empty()) return {};
    return result;
}

std::string choose_timeout_command() {
    std::string timeout_path = find_command_in_path("timeout");
    if (!timeout_path.empty()) return timeout_path;

    std::string gtimeout_path = find_command_in_path("gtimeout");
    if (!gtimeout_path.empty()) return gtimeout_path;

    return {};
}

// ---------------------------------------------------------------------------
// Parse a DIMACS solution (v-lines) from solver stdout
// ---------------------------------------------------------------------------

std::map<int, bool> parse_solution(const std::string& output) {
    std::map<int, bool> assignment;
    std::istringstream stream(output);
    std::string line;
    while (std::getline(stream, line)) {
        if (line.empty()) continue;
        if (line[0] == 'v' && (line.size() == 1 || line[1] == ' ')) {
            std::istringstream ls(line.substr(1));
            int lit;
            while (ls >> lit) {
                if (lit == 0) break;
                assignment[std::abs(lit)] = (lit > 0);
            }
        }
    }
    return assignment;
}

// ---------------------------------------------------------------------------
// Parse factorization metadata from _cnf_solution.txt
// ---------------------------------------------------------------------------

struct BitMapping {
    int bit_index = 0;
    bool is_pinned = false;
    int pinned_value = 0;
    int cnf_var = 0;
    bool negated = false;
};

struct FactorizationMeta {
    long long p = 0, q = 0, N = 0;
    std::vector<BitMapping> p_bits;
    std::vector<BitMapping> q_bits;
    bool valid = false;
};

FactorizationMeta load_factorization_meta(const std::string& path) {
    FactorizationMeta meta;
    std::ifstream f(path);
    if (!f.is_open()) return meta;

    std::string line;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == 'c') continue;
        std::istringstream ls(line);
        std::string tag;
        ls >> tag;
        if (tag == "p") {
            ls >> meta.p;
        } else if (tag == "q") {
            ls >> meta.q;
        } else if (tag == "N") {
            ls >> meta.N;
        } else if (tag == "p_bit" || tag == "q_bit") {
            BitMapping bm;
            std::string kind;
            ls >> bm.bit_index >> kind;
            if (kind == "pinned") {
                bm.is_pinned = true;
                ls >> bm.pinned_value;
            } else if (kind == "var") {
                bm.is_pinned = false;
                ls >> bm.cnf_var;
                std::string pol;
                ls >> pol;
                bm.negated = (pol == "neg");
            } else {
                continue;
            }

            if (tag == "p_bit") meta.p_bits.push_back(bm);
            else meta.q_bits.push_back(bm);
        }
    }

    meta.valid = (meta.N > 0 && !meta.p_bits.empty() && !meta.q_bits.empty());
    return meta;
}

long long reconstruct_factor(const std::vector<BitMapping>& bits,
                             const std::map<int, bool>& assignment) {
    long long val = 0;
    for (const auto& bm : bits) {
        bool bit = false;
        if (bm.is_pinned) {
            bit = (bm.pinned_value != 0);
        } else {
            auto it = assignment.find(bm.cnf_var);
            if (it == assignment.end()) return -1;
            bit = it->second;
            if (bm.negated) bit = !bit;
        }
        if (bit) val |= (1LL << bm.bit_index);
    }
    return val;
}

// ---------------------------------------------------------------------------
// Verify a solver's assignment against the CNF instance
// ---------------------------------------------------------------------------

bool verify_against_cnf(const std::string& cnf_path,
                        const std::map<int, bool>& assignment) {
    std::ifstream f(cnf_path);
    if (!f.is_open()) return false;

    std::string line;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == 'c' || line[0] == 'p') continue;

        std::istringstream ls(line);
        int lit;
        bool satisfied = false;
        while (ls >> lit) {
            if (lit == 0) break;
            int var = std::abs(lit);
            auto it = assignment.find(var);
            if (it != assignment.end()) {
                bool val = it->second;
                if (lit < 0) val = !val;
                if (val) {
                    satisfied = true;
                    break;
                }
            }
        }
        if (!satisfied) return false;
    }
    return true;
}

// ---------------------------------------------------------------------------
// Run a solver as a subprocess with timeout
// ---------------------------------------------------------------------------

SolverRun run_solver(const std::string& solver_path,
                     const std::string& solver_name,
                     const std::string& cnf_path,
                     int timeout_sec,
                     const std::string& planted_path,
                     const std::string& timeout_cmd) {
    SolverRun run;
    run.solver_name = solver_name;
    run.instance = fs::path(cnf_path).filename().string();

    std::string solver_flags;
    if (solver_name == "glucose") {
        solver_flags = " -model";
    }

    std::string tmpfile = "/tmp/sat_bench_" + solver_name + "_" +
                          std::to_string(getpid()) + "_" +
                          std::to_string(std::chrono::high_resolution_clock::now().time_since_epoch().count()) +
                          ".out";

    std::string cmd;
    if (!timeout_cmd.empty()) {
        cmd = shell_quote(timeout_cmd) + " " + std::to_string(timeout_sec) + " " +
              shell_quote(solver_path) + solver_flags + " " +
              shell_quote(cnf_path) + " > " + shell_quote(tmpfile) + " 2>/dev/null";
    } else {
        cmd = shell_quote(solver_path) + solver_flags + " " + shell_quote(cnf_path) +
              " > " + shell_quote(tmpfile) + " 2>/dev/null";
    }

    auto t0 = std::chrono::high_resolution_clock::now();
    int ret = std::system(cmd.c_str());
    auto t1 = std::chrono::high_resolution_clock::now();
    run.elapsed_sec = std::chrono::duration<double>(t1 - t0).count();

    std::string contents = read_text_file(tmpfile);

    if (ret == -1) {
        run.result = Result::ERROR;
        std::remove(tmpfile.c_str());
        return run;
    }

    int exit_code = -1;
    if (WIFEXITED(ret)) {
        exit_code = WEXITSTATUS(ret);
    } else {
        run.result = Result::ERROR;
        std::remove(tmpfile.c_str());
        return run;
    }

    if (exit_code == 10) {
        run.result = Result::SAT;
    } else if (exit_code == 20) {
        run.result = Result::UNSAT;
    } else if (exit_code == 124) {
        run.result = Result::TIMEOUT;
    } else {
        if (contents.find("s SATISFIABLE") != std::string::npos) {
            run.result = Result::SAT;
        } else if (contents.find("s UNSATISFIABLE") != std::string::npos) {
            run.result = Result::UNSAT;
        } else if (exit_code == 0) {
            if (contents.find("\nv ") != std::string::npos ||
                contents.rfind("v ", 0) == 0) {
                run.result = Result::SAT;
            } else {
                run.result = Result::ERROR;
            }
        } else {
            run.result = Result::ERROR;
        }
    }

    if (run.result == Result::SAT) {
        auto assignment = parse_solution(contents);
        if (!assignment.empty()) {
            run.verification_attempted = true;
            run.solution_verified = verify_against_cnf(cnf_path, assignment);

            if (!planted_path.empty()) {
                auto meta = load_factorization_meta(planted_path);
                if (meta.valid) {
                    run.factorization_attempted = true;
                    run.recovered_p = reconstruct_factor(meta.p_bits, assignment);
                    run.recovered_q = reconstruct_factor(meta.q_bits, assignment);
                    if (run.recovered_p > 0 && run.recovered_q > 0) {
                        run.factorization_verified =
                            (run.recovered_p * run.recovered_q == meta.N);
                    }
                }
            }
        }
    }

    std::remove(tmpfile.c_str());
    return run;
}

// ---------------------------------------------------------------------------
// Discover solver executables
// ---------------------------------------------------------------------------

struct SolverInfo {
    std::string name;
    std::string path;
};

std::vector<SolverInfo> discover_solvers(const fs::path& exe_dir,
                                         const std::string& filter) {
    std::vector<std::pair<std::string, std::vector<std::string>>> candidates = {
        {"kissat",  {"kissat"}},
        {"cadical", {"cadical"}},
        {"glucose", {"glucose", "glucose-syrup", "glucose_static"}},
    };

    std::vector<SolverInfo> found;
    for (const auto& entry : candidates) {
        const auto& name = entry.first;
        const auto& binaries = entry.second;
        if (!filter.empty() && filter != name) continue;

        bool located = false;
        for (const auto& bin : binaries) {
            fs::path p = exe_dir / bin;
            if (fs::exists(p) && fs::is_regular_file(p)) {
                found.push_back({name, p.string()});
                located = true;
                break;
            }
        }
        if (located) continue;

        for (const auto& bin : binaries) {
            std::string result = find_command_in_path(bin);
            if (!result.empty() && fs::exists(result)) {
                found.push_back({name, result});
                located = true;
                break;
            }
        }
    }
    return found;
}

// ---------------------------------------------------------------------------
// Discover .cnf instance files
// ---------------------------------------------------------------------------

std::vector<std::string> discover_instances(const fs::path& dir,
                                            const std::string& single_file) {
    if (!single_file.empty()) {
        if (fs::exists(single_file)) {
            return {single_file};
        }
        std::cerr << "Error: file not found: " << single_file << "\n";
        return {};
    }

    std::vector<std::string> instances;
    for (const auto& entry : fs::directory_iterator(dir)) {
        if (entry.is_regular_file() && entry.path().extension() == ".cnf") {
            instances.push_back(entry.path().string());
        }
    }
    std::sort(instances.begin(), instances.end());
    return instances;
}

// ---------------------------------------------------------------------------
// Parse instance header for metadata
// ---------------------------------------------------------------------------

struct InstanceInfo {
    std::string filename;
    int num_vars = 0;
    int num_clauses = 0;
    std::string description;
};

InstanceInfo parse_instance_header(const std::string& path) {
    InstanceInfo info;
    info.filename = fs::path(path).filename().string();

    std::ifstream f(path);
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        if (line[0] == 'c') {
            if (info.description.empty() && line.size() >= 2) {
                info.description = line.substr(2);
            }
        } else if (line[0] == 'p') {
            std::istringstream ls(line);
            std::string p, cnf;
            ls >> p >> cnf >> info.num_vars >> info.num_clauses;
            break;
        }
    }
    return info;
}

// ---------------------------------------------------------------------------
// Build result tag string (same logic as before)
// ---------------------------------------------------------------------------

std::string make_result_tag(const SolverRun& run) {
    std::string tag = result_str(run.result);
    if (run.factorization_attempted) {
        tag += run.factorization_verified ? " F" : " F!";
    } else if (run.verification_attempted) {
        tag += run.solution_verified ? " C" : " C!";
    }
    return tag;
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------

int main(int argc, char* argv[]) {
    int timeout_sec = 60000000;
    std::string solver_filter;
    std::string single_file;
    std::string tsv_path;
    bool tsv_explicit = false;

    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if ((arg == "--timeout" || arg == "-t") && i + 1 < argc) {
            timeout_sec = std::atoi(argv[++i]);
        } else if ((arg == "--solver" || arg == "-s") && i + 1 < argc) {
            solver_filter = argv[++i];
        } else if ((arg == "--file" || arg == "-f") && i + 1 < argc) {
            single_file = argv[++i];
        } else if ((arg == "--tsv" || arg == "-c") && i + 1 < argc) {
            tsv_path = argv[++i];
            tsv_explicit = true;
        } else if (arg == "--help" || arg == "-h") {
            std::cout << "Usage: " << argv[0]
                      << " [OPTIONS] [INSTANCE.cnf]\n\n"
                      << "  INSTANCE.cnf    Run on a single .cnf file (positional argument)\n"
                      << "  --timeout SEC   Per-solver timeout in seconds (default: 60)\n"
                      << "  --solver NAME   Run only this solver (kissat|cadical|glucose)\n"
                      << "  --file FILE     Run on a single .cnf file (same as positional)\n"
                      << "  --tsv FILE      Write results to TSV file (overrides auto-naming)\n\n"
                      << "When an instance file is given (positional or --file), the output\n"
                      << "TSV is automatically named res_<basename>.tsv in the instance's\n"
                      << "directory, unless --tsv is specified explicitly.\n";
            return 0;
        } else if (arg[0] != '-') {
            // Positional argument: treat as instance file
            single_file = arg;
        }
    }

    // Auto-derive output TSV path from instance filename if not explicit
    if (!tsv_explicit && !single_file.empty()) {
        fs::path inst(single_file);
        std::string stem = inst.stem().string();
        fs::path parent = inst.parent_path();
        if (parent.empty()) parent = ".";
        tsv_path = (parent / ("res_" + stem + ".tsv")).string();
    }

    fs::path exe_path;
    try {
        exe_path = fs::canonical(fs::path(argv[0]));
    } catch (...) {
        exe_path = fs::absolute(fs::path(argv[0]));
    }
    fs::path exe_dir = exe_path.parent_path();

    auto solvers = discover_solvers(exe_dir, solver_filter);
    if (solvers.empty()) {
        std::cerr << "No SAT solvers found. Place kissat, cadical, and/or glucose\n"
                  << "in the same directory as this executable, or ensure they are in PATH.\n";
        return 1;
    }

    std::string timeout_cmd = choose_timeout_command();

    std::cout << "=== SAT Benchmark ===\n";
    std::cout << "Timeout: " << timeout_sec << "s\n";
    if (!timeout_cmd.empty()) {
        std::cout << "Timeout command: " << timeout_cmd << "\n";
    } else {
        std::cout << "Timeout command: none found (running without external timeout enforcement)\n";
    }
    std::cout << "Solvers: ";
    for (size_t i = 0; i < solvers.size(); i++) {
        if (i > 0) std::cout << ", ";
        std::cout << solvers[i].name << " (" << solvers[i].path << ")";
    }
    std::cout << "\n";

    auto instances = discover_instances(exe_dir, single_file);
    if (instances.empty()) {
        std::cerr << "No .cnf files found in " << exe_dir << "\n";
        return 1;
    }
    std::cout << "Instances: " << instances.size() << "\n\n";

    // --- Console table header (same as before) ---
    std::printf("%-40s %8s %8s | ", "Instance", "Vars", "Clauses");
    for (const auto& s : solvers) {
        std::printf("%10s %8s %-8s| ", s.name.c_str(), "Time(s)", "Result");
    }
    std::printf("\n");
    int line_len = 42 + 8 + 8 + 3 + static_cast<int>(solvers.size()) * (10 + 8 + 8 + 2);
    std::printf("%s\n", std::string(line_len, '-').c_str());

    // --- Open TSV output if requested ---
    std::ofstream tsv;
    if (!tsv_path.empty()) {
        tsv.open(tsv_path);
        if (!tsv.is_open()) {
            std::cerr << "Warning: could not open TSV output file: " << tsv_path << "\n";
        } else {
            // Write TSV header
            tsv << "instance\tnum_vars\tnum_clauses";
            for (const auto& s : solvers) {
                tsv << "\t" << s.name << "_time\t" << s.name << "_result";
            }
            tsv << "\n";
        }
    }

    for (const auto& inst_path : instances) {
        auto info = parse_instance_header(inst_path);

        std::string planted_path;
        {
            std::string base = inst_path;
            if (base.size() >= 4 && base.substr(base.size() - 4) == ".cnf") {
                base = base.substr(0, base.size() - 4);
            }
            std::string candidate = base + "_cnf_solution.txt";
            if (fs::exists(candidate)) planted_path = candidate;
        }

        // Console output
        std::printf("%-40s %8d %8d | ",
                    info.filename.c_str(), info.num_vars, info.num_clauses);
        std::fflush(stdout);

        // Start TSV row
        std::ostringstream tsv_row;
        if (tsv.is_open()) {
            tsv_row << info.filename << "\t" << info.num_vars << "\t" << info.num_clauses;
        }

        for (const auto& solver : solvers) {
            auto run = run_solver(solver.path, solver.name, inst_path,
                                  timeout_sec, planted_path, timeout_cmd);

            std::string result_tag = make_result_tag(run);

            // Console output
            std::printf("%10s %8.3f %-8s| ",
                        solver.name.c_str(), run.elapsed_sec, result_tag.c_str());
            std::fflush(stdout);

            // TSV columns for this solver
            if (tsv.is_open()) {
                tsv_row << "\t" << run.elapsed_sec << "\t" << result_tag;
            }
        }

        std::printf("\n");

        // Flush TSV row
        if (tsv.is_open()) {
            tsv << tsv_row.str() << "\n";
        }
    }

    std::printf("\n");
    std::printf("Result key: SAT/UNSAT/TIMEOUT/ERROR\n");
    std::printf("  F  = factorization verified (solver's assignment reconstructs p*q = N)\n");
    std::printf("  F! = factorization FAILED\n");
    std::printf("  C  = clauses verified (assignment satisfies all CNF clauses)\n");
    std::printf("  C! = clause verification FAILED\n");

    if (tsv.is_open()) {
        std::cout << "\nTSV results written to: " << tsv_path << "\n";
    }

    return 0;
}
