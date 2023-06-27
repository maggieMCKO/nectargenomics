#!/bin/bash

source "0_setupPhyloAcc.sh"

RUN_num="Target_NeFr40"
logDIR="log2_${RUN_num}"

targets=$(find $logDIR | grep "11426205" | grep "out")

for target in ${targets};
do
  if grep -q "time used:" ${target}; then
	   :  # : do nothing; if found this line, task was finished
  else
    out=$(basename ${target})
    echo -e "not finished: ${out}"
  fi

done

find $logDIR -maxdepth 0 | grep ".out" | grep -r "Elapsed" > step2_NeFr40.Elapsedtime
# cut -d ':' -f 3 step2_NeFr40.Elapsedtime | cut -d ' ' -f 2 | sort -u


# bash 2b_run_PhyloAcc_check.sh &> step2_NeFr40.check
# cat step2_NeFr40.check | sort