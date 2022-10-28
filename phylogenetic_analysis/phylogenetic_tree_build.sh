#!/bin/bash
#SBATCH -n 4
#SBATCH -t 4-00:00:00
#SBATCH --mem 96Gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=
#SBATCH -o %A_tp_wgs.out

####################
module load gubbins
module load mafft
module load raxml
###################

input_name=


##generate multiple alignment file from genome fastas with references included
mafft --memsave --maxiterate 2 ${input_name}.fasta > ${input_name}.aln.fasta

##run gubbins to identify recombinant regions and filter them, builds ML tree from recombinant-filtered SNP only alignment
run_gubbins.py --tree-builder raxml --model GTRGAMMA --threads 4 --converge-method recombination --filter-percentage 95 --iterations 2  --outgroup CP007548.1,CP002374.1 ${input_name}.aln.fasta

##run raxml to generate 1000 rapid bootstraps on recombinant-filtered SNP only alignment intermediate output from gubbins 
raxml HPC -s ${input_name}.aln.filtered_polymorphic_sites.fasta -n 1000R-boots_on_${input_name} -m GTRGAMMA -T 4 -p 12367 -x 12456 -N 1000 -f a -o ${outgroup} -k

##Apply bootstraps to ML tree to annotate confidence at each node
raxml HPC -f b -t ${input_name}.aln.final_tree.tre -z RAxML_bootstrap.1000R-boots_on_${input_name} -n boot-supp_1000R_${input_name}_on-tree -m GTRGAMMA -T 4 -o ${outgroup} -k











