#!/bin/bash
#SBATCH --account=teaching
#SBATCH --reservation=teaching
##SBATCH --ntasks=4               # number of MPI processes
#SBATCH --cpus-per-task 3
#SBATCH --mem-per-cpu=4G      # memory; default unit is megabytes
#SBATCH --time=0-3:00           # time (DD-HH:MM)
#SBATCH --output=%x-%j.out

module load samtools

for files in *.sam
do
        time samtools view -bS ${files} -o ${files%.*}_no_mt.bam
done

