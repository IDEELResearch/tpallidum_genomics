# *Treponema pallidum* WGS pipeline
**University of North Carolina at Chapel Hill**

**Infectious Disease Epidemiology and Ecology Lab (IDEEL, ideelresearch.org) | PI: Parr**

**Date:** 11/17/2021

**Author(s):** Original code from Nick Brazeau and GATK Best Practices. Adapted by Wentao Chen (see Chen W et al. J Infect Dis 2021) to include filters and optimization developed by David Smajs' group (see Grillová et al. Front. Microbiol 2019). Optimized for UNC's Longleaf cluster by Fredrick Nindo.

## Background

**Objective:** Establish a strict pipeline for analysis of *T. pallidum* whole-genome sequence data alignment, variant calling, extraction of outer membrane protein genes.

_Notes:_
- The pipeline has 6 phases as noted below.

- Raw sequencing data from the current Agilent hybrid capture -> Illumina sequencing approach come in 3 reads: R1, R2 and R3. R2 is a short 10-mer molecular barcode that we are not currently using and can be discarded. R3 is then renamed R2 to maintain the conventional PE read format/nomenclature for use in the pipeline.


## **Phase 1: QC of raw reads**

This step requires that Trimommatic be installed and access to illumina adaptor sequences for trimming. Paired end reads are trimmed and can either be compressed or uncompressed for the next step in the pipeline. Trimming aims to remove poor quality reads and adaptors from the raw reads.

## **Phase 2: Host read and contaminant (non-T.pallidum sequence) filtration**

Trimmed reads are thereafter screened for traces of host genome sequences using bbmap at a threshold of 2 minimum hits against the combined file of human, animal, plant and  fungal ribo sequences (i.e. hg19_main_mask_ribo_animal_allplant_allfungus.fa.gz). This is followed by running repair.sh embedded in bbmap for reads that may have been tossed or broken. 

Additionally, bacterial contamination is screened for by performing a search of the host-free, clean sequences against the strainseeker database using seqtk. 

## **Phase 3: Alignment and quality assessment of the alignment to reference**
Clean PE reads are aligned to the REFERENCE TPA genome_using BWA. If the clade is known a priori, it is advisable to use SS14 refseq or Nichols refseq for samples/reads identified as SS14-like or Nichols-like, respectively, as there are slight differences in the synteny of genomes between these two clades. 

The output of the alignment is bam file that is processed futher in the subsequent steps for example:
  - The output bam is sorted  and some base filterng is additionally done using samtools.
  - Duplicates are sought and removed from the resulting bam using picard.
  - Indels that may have been introduced during alignment steps are removed using GATK based on the REF intervals used in the alignment.
  - Mismatches in the alignment are removed using bamutils in the first instance and picard in the second and third steps. 
  - If there are too many soft-clipped reads in the resultant bam, it is advisable to remove them using bamutils.
  - Similarly, hardclipped reads can be reomved using picard. 
  - 8 auxiliary shell scripts are provided to perform various postfiltering steps. For instance, they are useful in removing hard clipped reads, chimeric alignments, removing repetetively aligned reads using the mapping quality threshold of 10 (MAPQ10.sh) and retention of high quality mapping alignments using  (MAPQ40.sh) script.
  - If there are any singletons after filtering, they should be removed using 8.remove_singletons.sh script.
  - Before variant calling phase, it is useful to perform alignment quality ie bam QC checks using qualimap.

## **Phase 4: Variant calling and filtration**

In this phase, quality-filtered alignment files from previous steps are used to call, genotype, filter and select high quality variants that will in the end be used for generating consensus whole genome fasta sequences in a manner similar to reference-based genome assembly.

This is achieved through the following steps:

 - We use GATK to perform haptotype calling. At the time of writing this pipeline, GATK version 3 was employed and therefore it is advised to make appropriate modifications for calling certain GATK functions if newer versions are used. This step will generate a vcf file for each sample in the alignment.

 - To perform joint variant calling (ie genotype the haplotype), we manually combine each sample's raw genotype vcf file as in the two step commands below:
    - ls *vcf > vcfs.list
    - /nas/longleaf/apps/gatk/3.8-0/gatk -T CombineGVCFs -R /proj/ideel/resources/genomes/Tpallidum/SS14_CP004011.1.fasta --variant vcfs.list -o ss14_aligned.haplocaller.joint.raw.g.vcf

 - Once we have generated the combined raw gvcf file, we can now perform joint variant calling using the Genotype function in GATK. This will generate a joint raw vcf file. 

 - We subsequently filter variants based on set parameters to retain those that are reliable and of high quality. Some of the cut-offs and threshold levels are as shown in the command below:

    - /nas/longleaf/apps/gatk/3.8-0/gatk -T VariantFiltration -R {REF} -V {input} --filterExpression "QD < 2.0" --filterName "QD" --filterExpression "MQ < 40.0" --filterName "MQ" --filterExpression "FS > 60.0" --filterName "FS" --filterExpression "SOR > 3" --filterName "SORhigh" --filterExpression "MQRankSum <  -12.5" --filterName "MQRankSum" --filterExpression "ReadPosRankSum < -8.0 " --filterName "ReadPosRankSum" --filterExpression "DP < 3" --filterName "LowCoverage"

- We further select variants that meet the filtering criteria above - ie retain those that passed the quality filtering step in the previous step.

- This file is used as input for the next step that makes use of vcftools to generate a SNP only vcf, used to generate a consensus whole-genome fasta from vcf file.


## **Phase 5: Generation of whole-genome consensus fasta files**

Once you generate the SNP-only vcf, you must run the GetTPseq_fix_replaceDot.java (java script) to get the consensus fasta sequence. 

-This phase makes use of read depth information. Read depth can be obtained by running this command:

 samtools depth input.bam > input.bam.depth.txt
 
-Before running the java script, you need to set up the file PATH for the GetTPseq_fix_replaceDot.java file, for example:

 export PATH=/proj/ideel/jonbparrlab/users/fnindo/test10/all_runs/GetTPseq_fix_replaceDot.java:$PATH
 
-run (compile) the java script by running the command:

 javac GetTPseq_fix_replaceDot.java
 
-Generate the consensus fastas from the SNP-only vcf by running this command:

 java GetTPseq_fix_replaceDot


## **Phase 6: Extraction of OMPeome sequences of interest**

To make an OMPeome interval file that facilitates the extraction of OMPeome sequences from consensus whole genome fastas (generated in Phase 5) given their seqIDs in a text file and intervals_file, run the following command:

while IFS= read -r line1; do while IFS= read -r line2; do echo "$line1$line2"; done <test_intervals.txt; done <nichols_clade1.txt
