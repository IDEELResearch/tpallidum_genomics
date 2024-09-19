## *Treponema pallidum* Phylogenetic analysis

#### 1- Generate multiple alignment file from consensus genomes.
We used mafft to create a multiple alignment file from the variant called consensus genome output of the illumina_pipeline. Multiple reference genomes were included in this step (Nichols, SS14, Mexico_A, as well as SamoaD(TPE) and BosniaA(TEN) for an outgroup).


#### 2- Identification of recombinant regions, and masking
Gubbins was used to identify recombinant regions of the genome. The multiple alignment file was filtered down to a SNP-only, recombinant masked, multiple alignment file. 

#### 3- Rapid bootstrap analysis
The filtered SNP alignment was run through raxML to generate 1000 rapid bootstrapped trees.

#### 4- Applying bootstrap support to best ML tree.
RaxML was used to apply bootstrap support data to best ML tree. 