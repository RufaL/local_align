#!/bin/bash
#SBATCH --job-name=570_job
#SBATCH --mail-user=rufal@umich.edu
#SBATCH --mail-type=END
#SBATCH --cpus=per-task=1
#SBATCH --nodes=1
#SBATCH --ntacks-per-node=1
#SBATCH --mem-per-cpu=1000m
#SBATCH --time=10:00
#SBATCH --account=eecs570w20_class
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --output=/home/rufal/genacc/local_algin/swalign-cu.log

# Load CUDA
module load cuda/10.1.243
nvcc --version

# Compile
cd ~/genacc/local_align
nvcc swalign.cu -o swout_cu

# Run
./swout_cu
