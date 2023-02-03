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


bamutils=/nas/longleaf/home/hennelly/ngsutils/bin/bamutils
#ngsutils=



############################################################################

# Things to change 

############################################################################

module load python/2.7.12


#array=($(ls $sample_list | cut -d'_' -f1))
array=($(cat $sample_list))

mkdir -p ${pipeline_dir}/postfilter
mkdir -p ${pipeline_dir}/depth

for i in "${array[@]}"; do


#mm1
${bamutils} filter ${pipeline_dir}/aln/${i}.realn.bam ${pipeline_dir}/postfilter/${i}.mm1.bam -failed ${pipeline_dir}/postfilter/${i}.mm1.fail.txt -maximum_mismatch_ratio 5

#mm2
${bamutils} filter ${pipeline_dir}/postfilter/${i}.mm1.bam ${pipeline_dir}/postfilter/${i}.mm2.bam -failed ${pipeline_dir}/postfilter/${i}.mm2.fail.txt -maximum_mismatch_ratio 5

#mm3
cat ${pipeline_dir}/postfilter/${i}.mm1.fail.txt ${pipeline_dir}/postfilter/${i}.mm2.fail.txt > ${pipeline_dir}/postfilter/${i}.mm.fail.txt

##proceed to next few rules: scripts 1 and then 2 

#########softclipping
#${bamutils} removeclipping ${pipeline_dir}/postfilter/${i}.mm2.bam.short.bam ${pipeline_dir}/postfilter/${i}.scf.bam.tmp



done







echo "Read through script."
