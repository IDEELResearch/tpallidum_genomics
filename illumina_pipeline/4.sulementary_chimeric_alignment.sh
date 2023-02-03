#!/bin/bash
#SBATCH -n 1
#SBATCH -t 11-00:00:00
#SBATCH --mem 4800
#SBATCH --mail-type=ALL
#SBATCH --mail-user=fredrick_nindo@med.unc.edu
#SBATCH -o %A_tp_wgs.out
module load samtools
for sample in *.scf.bam; do samtools view -@ 24 -h -f 2048 -b $sample > ${sample}.sup.filtOut.bam & samtools view -@ 24 -h -F 2048 -b $sample > ${sample}.sup.bam;done
