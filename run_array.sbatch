#!/bin/bash

#SBATCH --job-name=ps_replication_array
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --array=1-2
#SBATCH --cpus-per-task=30
#SBATCH --time=168:00:00
#SBATCH --mem=10GB
#SBATCH --output=ps_rep_%A_%a.out
#SBATCH --error=ps_rep_%A_%a.err

module purge
module load python/intel/3.8.6
python parallel_driver.py batch_$SLURM_ARRAY_TASK_ID.csv EXT_data
