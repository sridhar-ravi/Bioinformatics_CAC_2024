#!/bin/bash
##SBATCH --account=teaching # not required for Frontenac 2.0
##SBATCH --reservation=teaching # not required for Frontenac 2.0
##SBATCH --ntasks=4               # number of MPI processes
#SBATCH --cpus-per-task 24
#SBATCH --mem-per-cpu=4G      # memory; default unit is megabytes
#SBATCH --time=0-3:00           # time (DD-HH:MM)
#SBATCH --output=%x-%j.out

ml samtools

find . -name "*.sam" | parallel -j 4 "time samtools view -@ 6 -bS {} | samtools sort -o {.}_mt_sorted.bam"


