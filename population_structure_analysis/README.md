## *Treponema pallidum* population structure analysis

#### 1- Read alignment
For *T. pallidum* population structure analysis we mapped all the sequences to a modified Nichols reference genome (CP004010.2) masked for *arp*, *TP0470* and *tpr* genes (see variant_calling.sh).

#### 2- PCA analysis
For the PCA analysis we filtered the VCF file to only include SNVs that were biallelic in *T pallidum* population and have minor allele frequency (MAF) greater than 0.05. We also excluded SNVs with excessive missingness (more than 0.2 missing per SNVs). Finally, we converted the VCF file to EIGENSTRAT format.

#### 3- FastBAPS analysis
For the fastbaps we first constructed a maximum likelihood phylogeny using IQ-tree. Then BAPS were determined with fastbaps. 
