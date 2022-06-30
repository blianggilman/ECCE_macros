#!/bin/bash
#SBATCH --image=docker:reynier0611/rhic_sl7_ext
#SBATCH --qos=shared
#SBATCH --constraint=haswell
#SBATCH --time=1:00:00
#SBATCH --array=0-199
shifter ./shifter.sh $SLURM_ARRAY_JOB_ID $SLURM_ARRAY_TASK_ID 10000