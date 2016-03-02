#!/bin/bash

#SBATCH --job-name=EkmanBL
#SBATCH --output=EkmanBL-%j.out
#SBATCH --error=EkmanBL-%j.err
#SBATCH -p normal

#SBATCH --time=48:00:00
#SBATCH --qos=

#SBATCH --nodes=4
#SBATCH --ntasks-per-node=16

# Load required modules
module load openmpi/1.8.7/intel
cd /home/aditya90/Codes/PadeOps/build/problems/incompressible/
mpirun -np 64 ./EkmanBL EkmanBL_files/input_Ekman.dat
