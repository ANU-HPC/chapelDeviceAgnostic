#!/bin/bash

#PBS -P vp91
#PBS -q gpuvolta
#PBS -l ncpus=12
#PBS -l ngpus=1
#PBS -l mem=100GB
#PBS -l walltime=00:02:00
#PBS -l wd

module load cuda/12.3.2

./heat_gpu -sn=8000 -snsteps 20
