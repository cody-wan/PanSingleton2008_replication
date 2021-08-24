#!/bin/bash
# 
#SBATCH --job-name=ps_replication_batch_4
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=5GB
#SBATCH --time=168:00:00

module purge
module load python/intel/3.8.6
​python parallel_driver.py
