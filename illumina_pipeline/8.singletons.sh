#!/bin/bash
#SBATCH -n 1
#SBATCH -t 11-00:00:00
#SBATCH --mem 4800
#SBATCH --mail-type=ALL
#SBATCH --mail-user=fredrick_nindo@med.unc.edu
#SBATCH -o %A_tp_wgs.out
module load samtools

for sample in *.scf.bam.sup.bam.MAPQ10.bam.MAPQ40.bam; do samtools view -@ 24 $sample | cut -f 1 | sort -T /tmp --parallel=24 | uniq -u > ${sample}.singleAfterFilt.txt; done
