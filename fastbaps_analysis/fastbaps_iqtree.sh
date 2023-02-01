#!/bin/bash

#SBATCH -c 24
#SBATCH -t 10-0:00:00
#SBATCH --mem 96Gb
#SBATCH -o %A_iqtree.out


####################
module load iqtree/1.6.12
####################

# Running IQ-Tree to build Maximum Likelihood phylogenetic tree

iqtree -s U19_uniq-TPA_Public_nic_SNP_only.fasta -bb 1000 -nm 20000 -nt AUTO -safe
