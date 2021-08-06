# tpallidum_genomics
#Overview
#Develop a strict pipeline to analysis Treponema pallidum Whome genome sequences data.
# Author: Wentao Chen but adopted to run on unc longleaf cluster by Fredrick Nindo
# Date : 2021/12/12
# NOTE : Install the Tp analysis environment first. More details : https://github.com/opplatek/bacterial_genome_analysis
These pipeline has 4 major phases:

Phase 1 QC of raw reads

Phase 2 Host read and contaminant (non-T.pallidum sequence) filteration

Phase 3 Alignmnent and Quality assessment of the alignment to Reference

Phase 4 Variant calling and Filtration

Phase 5 Generation of whole Genome consensus fasta.
