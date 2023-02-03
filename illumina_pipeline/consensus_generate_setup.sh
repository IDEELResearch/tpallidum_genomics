#!/bin/bash

set -o errexit
set -o nounset
set -o pipefail
set -o xtrace


proj_dir=/pipeline
rawreads=/rawreads
sample_list=${proj_dir}/resources/sample_list.txt

analysis_run=ss14_pipeline

pipeline_dir=${proj_dir}/${analysis_run}
#sample_list=${pipeline_dir}/symlinks

pipeline_output=${proj_dir}/variants/HaploCaller_SNP_joint.pass.vcf



############################################################################

# Things to change 

############################################################################

module load python/2.7.12
module load samtools
module load vcftools

#array=($(ls $sample_list | cut -d'_' -f1))
array=($(cat $sample_list))

mkdir -p ${pipeline_dir}/postfilter
mkdir -p ${pipeline_dir}/depth

for i in "${array[@]}"; do


##generate samtools depth for java scripts at end of tp pipeline
samtools depth ${pipeline_dir}/postfilter/${i}.MAPQ40.singletons.final.bam > ${pipeline_dir}/depth/${i}.depth.txt


done

##remove indels from vcf output of pipeline
vcftools --vcf ${pipeline_output} --remove-indels --recode --recode-INFO-all --out SNPs_only

##normal list generate
cat $sample_list > ${pipeline_dir}/consensus_generate/sample_list_half.txt
##list formatted for consensus_generate java script
paste ${pipeline_dir}/consensus_generate/sample_list_half.txt ${pipeline_dir}/consensus_generate/sample_list_half.txt > ${pipeline_dir}/consensus_generate/samples.txt






echo "Read through script."
