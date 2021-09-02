#!/bin/bash
#SBATCH -n 1
#SBATCH -t 11-00:00:00
#SBATCH --mem 4800
#SBATCH --mail-type=ALL
#SBATCH --mail-user=fredrick_nindo@med.unc.edu
#SBATCH -o %A_tp_wgs.out
module load samtools
for sample in *sup.bam; do samtools view -@ 24 -h $sample | awk '$6 ~/H/{print} ' | samtools view -bh - > ${sample}.hcf.filtOut.bam & samtools view -@ 24 -h $sample | awk ' 6 !~ /H/{print} ' | samtools view -bh - > ${sample}.hcf.bam; done 
