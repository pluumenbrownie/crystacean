#!/bin/bash

#SBATCH --nodes 1

#SBATCH --time 96:00:00 

#SBATCH --ntasks 32
 
export OMP_NUM_THREADS=1
 
module load mpi
 
srun /mnt/scratch/Work/jcottom/software/cp2k/exe/local/cp2k.popt testlarge.in  # <<== use srun and there is no need for mpiexec or mpirun or any -n / -np options 
