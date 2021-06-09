#!/bin/bash
#SBATCH --partition=singleGPU
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=2-00:00:00
#SBATCH --job-name=simcov_test_w_model
#SBATCh --output=simcov.out

module load upcxx/2020.10.0-python3-3o75
module load cmake/3.18.4-2lmi
export UPCXX_NETWORK=ibv

#cd lungmodel
#rm -rf lung_model_data
#./lungmodel --dim 300 300 300 --levels 0 --scale 10 --output lung_model_data
#./lungmodel --dim 300 300 300 --levels 3 --scale 2000 --output lung_model_data
#./lungmodel --dim 800 800 800 --levels 10 --scale 2000 --output lung_model_data
#cd ..
#upcxx-run -N 1 -n 16 -- install/bin/simcov --config covid_default.config --progress -v
./lungmodel
