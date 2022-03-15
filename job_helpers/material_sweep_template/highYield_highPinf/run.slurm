#!/bin/bash
#
#SBATCH --job-name=LAD_tuning
#SBATCH --time=8:00:00
#SBATCH --output=job-%j.out
#SBATCH --error=job-%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mail-user=jrwest@stanford.edu
#SBATCH --mail-type=ALL
#SBATCH --partition=normal
##########################
source ~/.bashrc
setup_padeops
ulimit -s unlimited

run_directory=`pwd`
rm velocityDiscontinuity

cd /home/users/jrwest/Codes/PadeOps
rm -r temp_build_9
mkdir temp_build_9
cd temp_build_9
cmake ..
cd problems/solids
make velocityDiscontinuity
cp velocityDiscontinuity $run_directory 
cd $run_directory

for f in ./input_C*.dat 
do 
    echo "============================================"
    echo "============================================"
    echo "============ Starting new job: ============="
    echo "============================================"
    echo "============================================"
    mpiexec -n 1 ./velocityDiscontinuity $f
done

echo "Done"
