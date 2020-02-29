#!/bin/bash -l

# Standard output and error:
#SBATCH -o ./tjob.out.%j
#SBATCH -e ./tjob.err.%j

# Initial working directory:
#SBATCH -D ./

# Job Name:
#SBATCH -J gDisk

# Queue (Partition):



# Number of nodes and MPI tasks per node:
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4

# Enable Hyperthreading
#SBATCH --ntasks-per-core=2

# for OpenMP/MPI?:
#SBATCH --cpus-per-task=8

# Request size of main memory per node in units of MB:
#SBATCH --mem=10240


# Run the program
python3 make_IC.py
srun ../GIZMO my.params > prog.out
