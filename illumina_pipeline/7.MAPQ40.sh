#!/bin/bash
#SBATCH -n 1
#SBATCH -t 11-00:00:00
#SBATCH --mem 4800
#SBATCH --mail-type=ALL
#SBATCH --mail-user=fredrick_nindo@med.unc.edu
#SBATCH -o %A_tp_wgs.out
module load samtools
for sample in *.scf.bam.sup.bam.MAPQ10.bam; do samtools view -@ 24 -h $sample | awk -v var="40" '$5 < var || $1 ~/^@/' | samtools view -b - > ${sample}.MAPQ40.filtOut.bam & samtools view -@ 24 -h -bq 40 $sample > ${sample}.MAPQ40.bam;done
