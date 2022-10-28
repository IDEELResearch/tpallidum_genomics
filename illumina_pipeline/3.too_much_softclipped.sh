#!/bin/bash
#SBATCH -n 1
#SBATCH -t 11-00:00:00
#SBATCH --mem 4800
#SBATCH --mail-type=ALL
#SBATCH --mail-user=fredrick_nindo@med.unc.edu
#SBATCH -o %A_tp_wgs.out

module load samtools

for sample in *scf.bam.tmp;
do samtools view -@ 24 $sample | grep 'ZC:f:' | awk '{for (i=1;i<=NF;i++){if ($sample ~/ZC:f:/) {print $1, $sample}}}' | sed 's/ZC:f://' | awk -v MAX_SOFTCLIP="0.05" '$2 > MAX_SOFTCLIP {print $1}' > ${sample}.read_names_to_remove_highSoftClip.txt.tmp # Get only softclipped reads above the threshold
    cat ${sample}.read_names_to_remove_highSoftClip.txt.tmp | cut -f 1 | sort -T /tmp --parallel=24 | uniq > ${sample}.read_names_to_remove_highSoftClip.txt;done # Get unique read names with softclipping!
    #cat ${sample%.*}.read_names_to_remove_highSoftClip.txt | cut -f 1 | sort -T $OUTPUT_DIR/tmp --parallel=$THREADS | uniq -d > tmp; mv tmp ${sample%.*}.read_names_to_remove_highSoftClip.txt # Get unique read names with softclipping if whole pair failed the filtering!
