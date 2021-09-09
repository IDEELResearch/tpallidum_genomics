# tpallidum_genomics
#Overview
#Develop a strict pipeline to analysis Treponema pallidum Whome genome sequences data.
# Author: Wentao Chen but adopted to run on UNC longleaf cluster by Fredrick Nindo
# Date : 2020/12/12
These pipeline has 4 major phases:
Raw reads from sequencing platforms come in 3 strands, R1, R2 and R3. R2 is short oligo strand that should be discarded and R3 renamed R2 to maintain the PE format of reads before creating symbolic links for the reads to the analysis working directory.
**Phase 1 QC of raw reads**

This step requires that Trimommatic be installed and access to illumina adaptor sequences for trimming. Paired end reads are trimmed and can either be compressed or uncompressed for the next step in the pipeline. Trimming aims to remove poor quality reads and adaptors from the raw reads.

**Phase 2 Host read and contaminant (non-T.pallidum sequence) filteration**

-Trimmed reads are thereafter screened for traces of host genome sequences using bbmap at a threshold of 2 minimum hits against the combined file of human, animal, plant and  fungal ribo sequences (i.e. hg19_main_mask_ribo_animal_allplant_allfungus.fa.gz). This is followed by running repair.sh  embedded in bbmap for reads that may have been tossed or broken. 

-Additionally, bacterial contamination is removed by performing a search of the host-free clean sequences against the strainseeker database using seqtk. This step will generate clean paired end reads for subsequent steps of the pipeline.

**Phase 3 Alignmnent and Quality assessment of the alignment to Reference**
- Clean PE reads are aligned to the REFERENCE TPA genome_using BWA. If the clade is known a priori, it is advisable to use SS14 refseq and Nichols refseq for samples/reads identified as SS14-like or Nichols-like respectively as there are slight differences in the synteny of genomes between these two clades. 

- The output of the alignment is bam file that is processed futher in the subsequent steps for example:

  - the output bam is sorted  and some base filterng is additionally done using samtools.
  
  - duplicates are sought and removed from the resuling bam using picard.
  
  - Indels that may have been introduced during alignment steps are removed using GATK based on the REF intervals used in the alignment.
  
  - mismatches in the alignment are removed using bamutils in the first instance and picard in the second and third steps. 
  
  - if there are too much soft-clipped reads in the resultant bam, it is advisable to remove them using bamutils.

  - similarly hardclipped reads can be reomved using picard. 
  
  - 8 auxiliary shell scripts are provided to perform various postfiltering steps. for instance they are useeful in removing hard clipped reads, chimeric alignments, removing repetetively aligned reads usiing the mapping quality threshold of 10 (MAPQ10.sh) and retention of high quality mapping alignments at   MAPQ40.sh)

  - If there are any singletons after filtering, they should be removed using 8.remove_singletons.sh script.

  - Before variant calling phase it is useful to perform alignment quality ie bam QC checks using qualimap.

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

