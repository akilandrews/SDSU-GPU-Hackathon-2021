#!/bin/bash
#SBATCH --partition=bigmem-3TB
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=2-00:00:00
#SBATCH --job-name=simcov_test_w_model
#SBATCH --mail-user=akilandrews@unm.edu
#SBATCH --mail-type=END

module load upcxx/2020.10.0-python3-3o75
module load cmake/3.18.4-2lmi
export UPCXX_NETWORK=ibv
./lungmodel 25255 21031 43734 0 0 0
