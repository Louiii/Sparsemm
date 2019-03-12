#!/bin/bash
#SBATCH -p test.q
#SBATCH -n 1
module purge
module load likwid/4.1
module load slurm/current
mtrx1=$1
mtrx2=$2
mtrx3=$3
mtrx4=$4
mtrx5=$5
mtrx6=$6
likwid-perfctr -f -C 0 -g DATA -m ./sparsemm out.matrix ${mtrx1}.matrix ${mtrx2}.matrix ${mtrx3}.matrix ${mtrx4}.matrix ${mtrx5}.matrix ${mtrx6}.matrix
