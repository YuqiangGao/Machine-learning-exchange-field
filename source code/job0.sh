#!/bin/bash
#SBATCH -t 48:00:00
#SBATCH -p cms13
#SBATCH -N 1  --tasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --output=wf22N-%j.out
#SBATCH --error=wf22N-%j.err
#SBATCH --mail-type=all
#SBATCH --
module load  mkl/11.0.1 fftw/3.3.3 compiler/intel/13.1.117 gcc/10.1

export QUIP_ROOT=$(pwd)
export QUIP_ARCH=linux_x86_64_gfortran
make