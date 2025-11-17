#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 0-03:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem-per-cpu=150G

ct=$(awk -v  awkvar="${SLURM_ARRAY_TASK_ID}" 'NR==awkvar' ages.csv)
Rscript augur.R ${ct}
