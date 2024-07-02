#!/bin/bash

#PBS -P vp91
#PBS -q normal
#PBS -l ncpus=48
#PBS -l mem=100GB
#PBS -l walltime=00:02:00
#PBS -l wd

./pi_cpu
