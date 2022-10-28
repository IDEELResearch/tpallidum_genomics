#!/bin/bash
#SBATCH -n 12
#SBATCH -t 11-00:00:00
#SBATCH --mem 4800
#SBATCH --mail-type=ALL
#SBATCH --mail-user=fredrick_nindo@med.unc.edu
#SBATCH -o %A_tp_wgs.out
module load trimmomatic 
module load picard 
module load bwa
module load bbmap
module load seqtk
#module load gatk
module load qualimap/2.2.1
module load r/4.1.0
#source /nas/longleaf/home/nzabanyi/miniconda3/etc/profile.d/conda.sh
#conda activate ngsutils
snakemake -s tp_wgs_test1.py --cluster "sbatch -n24 -t 4-00:00:00 --mem 49152 -o Cluster_%A_tp_wgs.out" -j 120 
