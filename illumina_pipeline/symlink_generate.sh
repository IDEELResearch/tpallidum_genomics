#!/bin/bash

set -o errexit
set -o nounset
set -o pipefail
set -o xtrace

proj_dir=/pipeline
rawreads=
sample_list=${proj_dir}/resources/sample_list.txt

analysis_run=ss14_pipeline
#sequence_run_A=
#sequence_run_B=


read_dir=${rawreads}/run01/
pipeline_dir=${proj_dir}/${analysis_run}
symlink_dir=${proj_dir}/${analysis_run}/symlinks


#read_dir_A=${rawreads}/run01
#pipeline_dir_A=${proj_dir}/${sequence_run_A}
#symlink_dir_A=${proj_dir}/${sequence_run_A}/symlinks

#read_dir_B=${rawreads}/run02
#pipeline_dir_B=${proj_dir}/${sequence_run_B}
#symlink_dir_B=${proj_dir}/${sequence_run_B}/symlinks

############################################################################

# Things to change 

############################################################################

#module load perl
#module load r/4.1.0

#array=($(ls $sample_list | cut -d'.' -f1))
array=($(cat $sample_list))

mkdir -p $pipeline_dir
mkdir -p $symlink_dir



for i in "${array[@]}"; do


###for single sequenced direct to analysis libraries
ln -s ${read_dir}/${i}_*_R1_001.fastq.gz ${symlink_dir}/${i}_R1.fastq.gz
ln -s ${read_dir}/${i}_*_R3_001.fastq.gz ${symlink_dir}/${i}_R2.fastq.gz

##above but for alt read input
#ln -s ${read_dir}/${i}_1.fastq.gz ${symlink_dir}/${i}_R1.fastq.gz
#ln -s ${read_dir}/${i}_2.fastq.gz ${symlink_dir}/${i}_R2.fastq.gz



####for resequenced samples that need initial analysis before merging 
#ln -s ${read_dir_A}/${i}_*_R1_001.fastq.gz ${symlink_dir}/${i}A_R1.fastq.gz
#ln -s ${read_dir_A}/${i}_*_R3_001.fastq.gz ${symlink_dir}/${i}A_R2.fastq.gz

#ln -s ${read_dir_B}/sample_*-${i}_*_R1_001.fastq.gz ${symlink_dir}/${i}B_R1.fastq.gz
#ln -s ${read_dir_B}/sample_*-${i}_*_R3_001.fastq.gz ${symlink_dir}/${i}B_R2.fastq.gz

##once you have merged need to dump the proper sample name into a new symlinks folder for pipeline running purposes
#ln -s ${read_dir_A}/${i}_*_R1_001.fastq.gz ${symlink_dir}/${i}_R1.fastq.gz


#####



done




echo "Read through script."
