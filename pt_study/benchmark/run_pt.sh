#!/bin/bash
# run_pt.sh <cnf> <move> <timeout> <seed>
f=$1; move=$2; to=$3; seed=$4
m=$(awk '/^p cnf/{print $4}' $f)
if [ "$move" = uniform ]; then
  M=$(python3 -c "import math;print(max(8,math.ceil(48*math.sqrt($m/4698))))")
  /home/claude/pt/pt_cnf $f --move uniform --replicas $M --tmin 0.15 --tmax 1.2 --time $to --seed $seed | head -1
else
  /home/claude/pt/pt_cnf $f --move focused --replicas 16 --tmin 0.1 --tmax 0.6 --time $to --seed $seed | head -1
fi
