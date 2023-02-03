#!/bin/bash
#SBATCH -n 12
#SBATCH -t 2-00:00:00
#SBATCH --mem 64Gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=chris_hennelly@med.unc.edu
#SBATCH -o %A_tp_wgs.out
module load java
module load trimmomatic 
module load picard 
module load bwa
module load bbmap
module load seqtk
module load gatk/3.8-0
module load qualimap/2.2.1
module load r/4.1.0
#module load python/2.7.12
#source /nas/longleaf/home/nzabanyi/miniconda3/etc/profile.d/conda.sh
#source /nas/longleaf/home/hennelly/anaconda3/etc/profile.d/conda.sh
#conda activate ngsutils
snakemake -s tp_wgs_test1.py --cluster "sbatch -n24 -t 4-00:00:00 --mem 49152 -o Cluster_%A_tp_wgs.out" -j 120 
