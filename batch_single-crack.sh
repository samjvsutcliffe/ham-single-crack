#!/bin/bash

# Request resources:
#SBATCH -c 32 #entire node
#SBATCH --time=12:00:0  # 6 hours (hours:minutes:seconds)
#SBATCH --mem=32G      # 1 GB RAM
#SBATCH -p shared

module load gcc
module load aocl

echo "Running code"
rm output/*
rm outframes/*
#export GC_THREADS=32

sbcl --dynamic-space-size 32000  --disable-debugger --load "single-crack.lisp" --quit
#grep slurm-5714700.out.out -oP -e "Melt.*\K\d+\.\d+"
