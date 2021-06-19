#!/bin/bash
#SBATCH --partition=singleGPU
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=2-00:00:00
#SBATCH --job-name=simcov_test_w_model
#SBATCh --mail-user=akilandrews@unm.edu
#SBATCH --mail-type=END

module load upcxx/2020.10.0-python3-3o75
module load cmake/3.18.4-2lmi
export UPCXX_NETWORK=ibv
./lungmodel 300 300 300 0 0 0
