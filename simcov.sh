#!/bin/bash
#SBATCH --partition=dualGPU
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=2-00:00
#SBATCH --job-name=simcov_test_w_model
#SBATCH --mail-user=akilandrews@unm.edu
#SBATCH --mail-type=END

#./lungmodel 25255 21031 43734 0 0 0
./lungmodel 8000 8000 1 8627 7010 0