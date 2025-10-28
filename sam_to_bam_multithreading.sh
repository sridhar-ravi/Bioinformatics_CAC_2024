#!/bin/bash
##SBATCH --account=teaching #not required for Frontenac 2.0
##SBATCH --reservation=teaching # not required for Frontenac 2.0
#SBATCH --cpus-per-task 4
#SBATCH --mem-per-cpu=4G      # memory; default unit is megabytes
#SBATCH --time=0-3:00           # time (DD-HH:MM)
#SBATCH --output=%x-%j.out

ml samtools/1.18

for files in *.sam
do
time samtools view -@ 4 -b ${files} | samtools sort -o ${files%.*}_mt_sorted.bam
done

