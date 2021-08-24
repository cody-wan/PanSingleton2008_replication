#!/bin/bash

#SBATCH --job-name=test
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=1GB
#SBATCH --time=1:00:00

​
module purge
module load python/intel/3.8.6
​
python hpc_test.py