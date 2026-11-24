#!/bin/bash
# run every task with 2 workers, logging versions like the Slurm job does
cd /home/claude/rerun
n=0
while read -r c; do echo "$n|$c"; n=$((n+1)); done < commands.txt | \
xargs -P 2 -I{} bash -c 'i="${1%%|*}"; c="${1#*|}"; { echo "kissat_version: $(./kissat --version)"; echo "cadical_version: $(./cadical --version)"; echo "cpu: $(grep -m1 "model name" /proc/cpuinfo | cut -d: -f2)"; echo "command: $c"; eval "$c"; echo "exit_status: $?"; } > logs/out_local_$i.txt 2>&1' _ {}
echo ALLDONE > logs/ALLDONE
