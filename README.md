# tpallidum_genomics
#Overview
#Develop a strict pipeline to analysis Treponema pallidum Whome genome sequences data.
# Author: Wentao Chen but adopted to run on UNC longleaf cluster by Fredrick Nindo
# Date : 2020/12/12
These pipeline has 4 major phases:

**Phase 1 QC of raw reads**

This step requires that Trimommatic be installed and access to illumina adaptor sequences for trimming. Paired end reads are trimmed and can either be compressed or uncompressed for the next step in the pipeline.

**Phase 2 Host read and contaminant (non-T.pallidum sequence) filteration**

**Phase 3 Alignmnent and Quality assessment of the alignment to Reference**

**Phase 4 Variant calling and Filtration**

**Phase 5 Generation of whole Genome consensus fasta**

Once you generate the SNPs only vcf, you must run the GetTPseq_fix_replaceDot.java (java script) to get the consensus fasta sequence. 

-This phase makes use of read depth information. These can be obtained by running this command:

 samtools depth input.bam > input.bam.depth.txt
 
-Before running the java script, you need to set up the file PATH for GetTPseq_fix_replaceDot.java file for example :

 export PATH=/proj/ideel/jonbparrlab/users/fnindo/test10/all_runs/GetTPseq_fix_replaceDot.java:$PATH
 
-run (compile) the java script by running the command:

 javac GetTPseq_fix_replaceDot.java
 
-Generate the consensus fastas from you snp only vcf by running this command:

 java GetTPseq_fix_replaceDot

