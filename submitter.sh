#!/bin/bash
if [ $# -gt 0 ]; then
	echo "sbatch $1 | grep -P -o '\d+'"
	jobid=`sbatch $1 | grep -P -o '\d+'`
	echo $jobid

	if [ $# -gt 1 ]; then
		shift
		while [ "$1" != "" ]; do
			echo "sbatch --dependency=afterok:$jobid $1 | grep -P -o '\d+'" 
			jobid=`sbatch --dependency=afterok:$jobid $1 | grep -P -o '\d+'`
			echo $jobid
			shift
		done
	fi
fi
squeue -u $UID
