# Develop a strict pipeline to analysis Treponema pallidum Whome genome sequences data.
# Author: Wentao Chen but adopted to run on unc longleaf cluster by Fredrick Nindo
# Date : 2020/11/12
# NOTE : Install the Tp analysis environment first. More details : https://github.com/opplatek/bacterial_genome_analysis


####### Working Directory and Project Specifics ############
WORKDIR = '/pine/scr/n/z/nzabanyi/test14'
readWD = '/pine/scr/n/z/nzabanyi/test14/'
SAMPLES, = glob_wildcards(readWD + 'symlinks/{samp}_1.fastq.gz')
CLEANSAMPLES, = glob_wildcards(readWD + 'host_remove/{samp}.clean.fastq.gz')
BAMLIST, = glob_wildcards(readWD + 'variant/{samp}_HapoCaller.raw.g.vcf')


################   REFERENCE    ##############
REF = '/proj/ideel/resources/genomes/Tpallidum/SS14_CP004011.1.fasta'
#REF1= '/home/chenwt/Project/T.pallidum_WGS/test/postfilter/filt/final_filtered/consensus_seq/Tp_SS14_GCA_000410555_1.fasta'
GFF = '/proj/ideel/jonbparrlab/users/fnindo/test1/GCF_000604125.1_ASM60412v1_genomic.gff'
STRAIN_DB = '/pine/scr/n/z/nzabanyi/test14/ss_db_w32_4324/'
bamutils = '/nas/longleaf/home/nzabanyi/ngsutils/bin/bamutils'


######## Tools to Call #########
######## Always on #########
PICARD = '/nas/longleaf/apps/picard/2.2.4/picard-tools-2.2.4/picard.jar'
TMPDIR = '/proj/ideel/jonbparrlab/users/fnindo/test1/TMPDIR'
TRIMMOMATIC = '/nas/longleaf/apps/trimmomatic/0.36/Trimmomatic-0.36/trimmomatic-0.36.jar'
seeker = '/pine/scr/n/z/nzabanyi/test14/seeker.pl'
gatk = '/nas/longleaf/apps/gatk/3.8-0/gatk'


##########################################################################################

rule all:
     # input: expand('symlinks/{samp}_1.PAIREDtrimmomatictrimmed.fastq',samp = SAMPLES) 
     # input: expand('host_remove/{samp}_1.clean.fastq.gz',samp = SAMPLES)
     # input: expand('strainseeker/{samp}.seeker.txt',samp = CLEANSAMPLES)
     #input: expand('aln/{samp}.sorted.bam', samp =SAMPLES) 
     # input: expand('alignment_stats/{samp}.filt.flagstat', samp = SAMPLES)
      input: expand('aln/{samp}.realn.bam', samp = SAMPLES) 
      # input: expand('postfilter/{samp}.mm.fail.txt', samp = SAMPLES)
      # input: expand('postfilter/{samp}.mm.filtOut.bam', samp = SAMPLES)
      # input: expand('postfilter/{samp}.scf.bam.tmp', samp =SAMPLES)
      # input: expand('postfilter/{samp}.scf.filtOut.bam',samp =SAMPLES)
      # input: expand('postfilter/{samp}.scf.bam', samp =SAMPLES)
      # input: expand('postfilter/{samp}.singletons.filtOut.bam',samp =SAMPLES)
      # input: expand ('postfilter/filt/final_filtered/qc/{samp}/{samp}.qualimap.pdf', samp = SAMPLES)
      # input: expand('variants/{samp}_HaploCaller.raw.g.vcf', samp = SAMPLES)
      # input: expand ('variants/UNC_HaploCaller_joint.raw.vcf')
      # input: expand ('variants/HaploCaller_joint.qual.vcf')
      # input: expand ('variants/HaploCaller_SNP_joint.pass.vcf')
###############################################################################

## Generate consensus sequence based on VCF
## Output a new vcf file from the input vcf file that removes any indell sites
#   vcftools --vcf input_file.vcf --remove-indels --recode --recode-INFO-all --out SNPs_only
#   Once you get the SNPs_only.vcf, you could run the GetTPseq_fix_replaceDot.java(java script) to get the consensus sequence. Also, you need to set up the file PATH for GetTPseq_fix_replaceDot.java, then run "javac GetTPseq_fix_replaceDot.java" and "java GetTPseq_fix_replaceDot" step by step. NOTE: Prepare the reads depth information by running "samtools depth input.bam > input.bam.depth.txt"

rule select_variants :
        input:  'variants/HaploCaller_joint.qual.vcf'
        output:  'variants/HaploCaller_SNP_joint.pass.vcf'
        shell:  'java -jar {/nas/longleaf/apps/gatk/3.8-0/gatk} \
                -T SelectVariants \
                -R {REF} -V {input} -o {output} \
                -select "vc.isNotFiltered()" '

rule filter_variants :
        input:  'variants/UNCTP_HaploCaller_joint.raw.vcf'
        output:  'variants/HaploCaller_joint.qual.vcf'
        shell: '/nas/longleaf/apps/gatk/3.8-0/gatk -T VariantFiltration -R {REF} -V {input} --filterExpression "QD < 2.0" --filterName "QD" --filterExpression "MQ < 40.0" --filterName "MQ" --filterExpression "FS > 60.0" --filterName "FS" --filterExpression "SOR > 3" --filterName "SORhigh" --filterExpression "MQRankSum <  -12.5" --filterName "MQRankSum" --filterExpression "ReadPosRankSum < -8.0 " --filterName "ReadPosRankSum" --filterExpression "DP < 3" --filterName "LowCoverage" --logging_level ERROR -o {output}'

rule genotype_GVCFs:
        input: readWD + 'variants/UNC_Tp_variant_gvcfs.vcf'
        output: 'variants/UNC_Tp_HaploCaller_joint.raw.vcf'
        shell: '/nas/longleaf/apps/gatk/3.8-0/gatk -T GenotypeGVCFs -R {REF} --variant {input} -o {output}'

# generate a gvcfs.list manually.

rule haplotype_caller:
        input: '/proj/ideel/jonbparrlab/users/fnindo/test10/run6/postfilter/filt/final_filtered/{samp}.MAPQ40.singletons.final.bam',
        output: 'variants/{samp}_HaploCaller.raw.g.vcf'
        shell:'/nas/longleaf/apps/gatk/3.8-0/gatk -T HaplotypeCaller -R {REF} -I {input} -ploidy 1 -stand_call_conf 30 --emitRefConfidence GVCF -o {output} --unsafe'
  
##### QC #######
rule bam_qc:
        input:'postfilter/filt/final_filtered/{samp}.MAPQ40.singletons.final.bam'
        output:'postfilter/filt/final_filtered/qc/{samp}/{samp}.qualimap.pdf','postfilter/filt/final_filtered/qc/{samp}/'
        shell:'qualimap bamqc -bam {input} -nt 24 -c -outdir {output[1]} -outfile {output[0]} -outformat PDF'


## NOTE: You could run get_alignment_consensusSeq.sh for consensus sequence generation.

##########################
### Post-filter set  #####
## David's paper set the parameters: 
MAX_PERC_OF_MM='0.05' # Maximum percentage of mismatches compared to the read length - bad if we have to error-prone reads but helps to remove false-positives
MAX_NUMBER_OF_MM='5' # Maximum number of mismatches
MIN_LENGTH_MAPPED='35' # Remove mappings that mapped with too few bases (remove excessive soft-clipping)
MAX_SOFTCLIP='0.05' # Maximal percentage of the reads allowed to be soft-clipped
MAX_HARDCLIP='0.00' # Maximal percentage of the reads allowed to be hard-clipped
MAPQ='10' # Minimal MAPQ for (probably) repetitive regions
MAPQ_FINAL='40' # Minimal MAPQ for final results

###########################################################################
########### Run 8.Singletons.sh  NOTE: Only works for PE sequencing #######
###########################################################################

rule singletons_remove:
        input: 'postfilter/{samp}.scf.bam.sup.bam.MAPQ10.bam.MAPQ40.bam.singleAfterFilt.txt','postfilter/{samp}.scf.bam.sup.bam.MAPQ10.bam.MAPQ40.bam'
        output: 'postfilter/{samp}.singletons.filtOut.bam','postfilter/{samp}.MAPQ40.singletons.final.bam'
        shell: 'picard FilterSamReads I={input[1]} O={output[0]} READ_LIST_FILE={input[0]} FILTER=includeReadList & picard FilterSamReads I={input[1]} O={output[1]} READ_LIST_FILE={input[0]} FILTER=excludeReadList'

## Run 8.Singletons.sh first, then go to next step. ##########

##############################################################
##### Run 7. MAPQ40.sh Only high quality alignments###########
##############################################################

##############################################################
#### Run 6. MAPQ10.sh very often repetitive alignments########
##############################################################

##############################################################
################# Run 5. Hardclipped.sh ######################
##############################################################

##############################################################
###### Run 4. supplementary_chimeric alignment.sh ############
##############################################################

##############################################################
########### 3. Remove too much softclipped ###################
##############################################################

rule remove_too_much_softclipped3:
        input: 'postfilter/{samp}.mm2.bam.short.bam','postfilter/{samp}.scf.bam.tmp.read_names_to_remove_highSoftClip.txt'
        output: 'postfilter/{samp}.scf.bam'
        shell:'picard FilterSamReads I={input[0]} O={output} READ_LIST_FILE={input[1]} FILTER=excludeReadList'

rule remove_too_much_softclipped2:
        input: 'postfilter/{samp}.mm2.bam.short.bam','postfilter/{samp}.scf.bam.tmp.read_names_to_remove_highSoftClip.txt'
        output: 'postfilter/{samp}.scf.filtOut.bam'
        shell:'picard FilterSamReads I={input[0]} O={output} READ_LIST_FILE={input[1]} FILTER=includeReadList'
## run 3.too_much_oftcliped.sh before you go to next step.
rule remove_too_much_softclipped1:
        input:'postfilter/{samp}.mm2.bam.short.bam'
        output:'postfilter/{samp}.scf.bam.tmp',
        shell:'sh {bamutils} removeclipping {input} {output}'

##############################################################################
### Before run the step ,you should run 'sh 2.remove_short_alignment.sh' first.
##############################################################################

#############################################################
############## 1. remove too many mismatches ################
#############################################################
rule remove_too_many_mismatches4:
        input:'postfilter/{samp}.mm.fail.txt.tmp','aln/{samp}.realn.bam'
        output:'postfilter/{samp}.mm.filtOut.bam'
        shell:'picard FilterSamReads I={input[1]} O={output} READ_LIST_FILE={input[0]} FILTER=includeReadList'

# Run 1.remove_too_much_mismatches.sh, then go ahead. NOTE: Do all the *.sh in postfilter folder before process nest step


rule remove_too_many_mismatches3:
        input: 'postfilter/{samp}.mm1.fail.txt','postfilter/{samp}.mm2.fail.txt'
        output:'postfilter/{samp}.mm.fail.txt'
        shell: 'cat {input[0]} {input[1]} > {output}'
rule remove_too_many_mismatches2:
        input: 'postfilter/{samp}.mm1.bam'
        output:'postfilter/{samp}.mm2.bam','postfilter/{samp}.mm2.fail.txt'
        shell: 'sh {bamutils} filter {input} {output[0]} -failed {output[1]} -mismatch {MAX_NUMBER_OF_MM}'
rule remove_too_many_mismatches1:
        input:'aln/{samp}.realn.bam'
        output: 'postfilter/{samp}.mm1.bam','postfilter/{samp}.mm1.fail.txt'
        shell: 'sh {bamutils} filter {input} {output[0]} -failed {output[1]} -maximum_mismatch_ratio {MAX_PERC_OF_MM}'

############################
######## Alignment #########
############################

rule realn_indels:
	input: bam = 'aln/{samp}.matefixed.bam', chrs = 'intervals/Tp_SS14.intervals', targets = 'aln/{samp}.realigner.intervals', 
	output: 'aln/{samp}.realn.bam'
	shell: '/nas/longleaf/apps/gatk/3.8-0/gatk -T IndelRealigner \
		-R {REF} -I {input.bam} \
		-L {input.chrs} -targetIntervals {input.targets} \
		-o {output}' 

rule find_indels:
	input: bam = 'aln/{samp}.matefixed.bam', index = 'aln/{samp}.matefixed.bam.bai', chrs = 'intervals/Tp_SS14.intervals'
	output: 'aln/{samp}.realigner.intervals'
	shell: '/nas/longleaf/apps/gatk/3.8-0/gatk -T RealignerTargetCreator \
		-R {REF} -I {input.bam} \
		-L {input.chrs} \
                -o {output}'

rule index_dedups: 
	input: 'aln/{samp}.matefixed.bam'
	output: 'aln/{samp}.matefixed.bam.bai'
	shell: 'java -jar {PICARD} BuildBamIndex INPUT={input} OUTPUT={output} TMP_DIR={TMPDIR}'


rule make_matefixed:
        input: 'aln/{samp}.dedup.bam'
	output: 'aln/{samp}.matefixed.bam'
	shell: 'java -jar {PICARD} FixMateInformation INPUT={input} OUTPUT={output} TMP_DIR={TMPDIR}'
	
	
rule mark_dups:
	input: 'aln/{samp}.filt.bam'
	output:'aln/{samp}.dedup.bam','aln/{samp}.dedup.metrics'
	shell: 'java -jar {PICARD} MarkDuplicates \
		I={input} O={output[0]} \
		METRICS_FILE={output[1]} \
		TMP_DIR={TMPDIR} REMOVE_DUPLICATES=true \
		MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=100'
    
rule base_filtering2:
        input: 'aln/{samp}.filt.bam'
        output: 'alignment_stats/{samp}.filt.flagstat'
        shell: 'samtools index -@ 20 {input} | samtools flagstat {input} > {output}'

rule base_filtering1:
        input:'aln/{samp}.sorted.bam',
        output:'alignment_stats/{samp}.flagstat','aln/{samp}.filt.bam'
        shell: 'samtools index -@ 20 {input} | samtools flagstat {input} > {output[0]} & samtools view -@ 20 -h -F 12 -f 2 -F 256 -b {input} | samtools sort -n -@ 20 - | samtools fixmate -O bam - - | samtools sort -@ 20 - > {output[1]}'

rule sort_bam:
	input: 'aln/{samp}.raw.bam'
	output: 'aln/{samp}.sorted.bam'
	shell: 'samtools sort -@ 12 {input} -o {output} '

rule fastq_to_bam:
	input: 'host_remove/{samp}_1.clean.fastq.gz', 'host_remove/{samp}_2.clean.fastq.gz'
	output: 'aln/{samp}.raw.bam'
	shell: 'bwa mem {REF} {input[0]} {input[1]} -R "@RG\\tID:bwa\\tPL:illumina\\tLB:{wildcards.samp}_lib\\tSM:{wildcards.samp[0]}{wildcards.samp[1]}{wildcards.samp[2]}{wildcards.samp[3]}{wildcards.samp[4]}{wildcards.samp[5]}{wildcards.samp[6]}{wildcards.samp[7]}{wildcards.samp[8]}{wildcards.samp[9]}" | samtools view -Sb - > {output}'

rule bacteria contamination scan_2:
        input: 'strainseeker/{samp}.sub'
        output: 'strainseeker/{samp}.seeker.txt'
        shell: 'perl {seeker} -i {input} -d {STRAIN_DB} -o {output}'

rule bacteria contamination scan_1:
        input: 'host_remove/{samp}.clean.fastq.gz'
        output: 'strainseeker/{samp}.sub'
        shell: 'seqtk sample -s100 {input} 1000000 > {output}'


# Number of raw reads
#for sample in symlinks/*R1*gz
#do
#    echo "Raw read count for file echo $(basename $sample)"
#    echo $(zcat $sample | wc -l)/4 | bc
#done

# Number of the preprocessed reads
#for sample in symlinks/*R1*PAIREDtrimmomatictrimmed.fastq.gz
#do
#    echo "Preprocessed read count for file echo $(basename $sample)"
#    echo $(zcat $sample | wc -l)/4 | bc
#done

# Number of preprocessed and host-cleaned reads
#for sample in preprocessed/*R1*clean.fastq.gz
#do
#    echo "Cleaned read count for file echo $(basename $sample)"
#    echo $(zcat $sample | wc -l)/4 | bc
#done


##NOTE: make a fasqtc and multiqc again manually.
## shell: fastqc --threads 20 --outdir clean_fastqc preprocessed/*.clean.fastq.gz
## multiqc --outdir clean_fastqc clean_fastqc

rule get_clean_fastq:
    input: 'preprocessed/{samp}.clean.fastq.gz'
    output: 'host_remove/{samp}_1.clean.fastq.gz','host_remove/{samp}_2.clean.fastq.gz'
    shell: '/nas/longleaf/apps/bbmap/38.82/bbmap/repair.sh tossbrokenreads  in={input} out1={output[0]} out2={output[1]}'

rule host_genome_removal:
    input: 'symlinks/{samp}_1.PAIREDtrimmomatictrimmed.fastq','symlinks/{samp}_2.PAIREDtrimmomatictrimmed.fastq'
    output: 'preprocessed/{samp}.clean.fastq.gz','preprocessed/{samp}.dirty.fastq.gz'
    shell:'/nas/longleaf/apps/bbmap/38.82/bbmap/bbmap.sh threads=12 -Xmx25g minid=0.95 maxindel=3 bandwidthratio=0.16 bandwidth=12 quickmatch fast minhits=2 path=/proj/ideel/jonbparrlab/users/fnindo/test5/ unpigz pigz in={input[0]} in2={input[1]} outu={output[0]} outm={output[1]}'

rule trim_illumina_Adaptors_fastqs:
    input: readWD + 'symlinks/{samp}_1.fastq.gz', readWD + 'symlinks/{samp}_2.fastq.gz'
    output: 'symlinks/{samp}_1.PAIREDtrimmomatictrimmed.fastq', 'symlinks/{samp}_1.UNPAIREDtrimmomatictrimmed.fastq', 'symlinks/{samp}_2.PAIREDtrimmomatictrimmed.fastq', 'symlinks/{samp}_2.UNPAIREDtrimmomatictrimmed.fastq',  
    shell: 'java -jar /nas/longleaf/apps/trimmomatic/0.36/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 24 -trimlog symlinks/trim_log.txt {input[0]} {input[1]} {output[0]} {output[1]} {output[2]} {output[3]} ILLUMINACLIP:/nas/longleaf/apps/trimmomatic/0.36/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10:8:TRUE LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36'
     # The TRUE at the end keeps the paired end reads in R2
     # Want to align the PAIRED trimmed
     #adapterremoval, https://github.com/MikkelSchubert/adapterremoval

