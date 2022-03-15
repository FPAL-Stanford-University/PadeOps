for d in */ 
do 
    cd "$d"
    sbatch run.slurm	
    cd ../
done
