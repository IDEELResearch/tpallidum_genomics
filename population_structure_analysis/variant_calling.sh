#!/bin/bash
#SBATCH -J u19_variant_calling
#SBATCH -p general
#SBATCH -n 8 
#SBATCH -t 8:00:00 
#SBATCH --mem 64g  
#SBATCH -o u19_vcf.%A.%a.out # Output file name
#SBATCH -e u19_vcf.%A.%a.err # Error file name
#SBATCH --array 0-1500%20

module load trimmomatic/0.36 bbmap/38.82 seqtk bamutils r/4.1.0 bwa/0.7.17 samtools/1.8 gatk/3.8.1 picard/2.2.4

DIR="/PATH/TO/THE/WORKING/DIRECTORY"

mkdir -p ${DIR}/trim
mkdir -p ${DIR}/host_clean
mkdir -p ${DIR}/aln
mkdir -p ${DIR}/tmp
mkdir -p ${DIR}/postfilteri
mkdir -p ${DIR}/U19_bam
mkdir -p ${DIR}/U19_gvcf

trimDIR="${DIR}/trim"
cleanDIR="${DIR}/host_clean"
cleanDIR="${DIR}/host_clean"
alignDIR="${DIR}/aln"
alignTMP="${DIR}/tmp"
filterDIR="${DIR}/postfilter"
alignDIR="${DIR}/aln"
bamDIR="${DIR}/U19_bam"
gvcfDIR="${DIR}/U19_gvcf"

SAMPLE=($(cat ${DIR}/u19_wgs.sampleList.txt))
REF="/TPallidum/Pallidum_Nichols/TPal_Nichols_CP004010.2_masked.fasta"

## trimming reads with trimmomatic
 
trimmomatic PE -threads 8 -trimlog ${trimDIR}/${SAMPLE[${SLURM_ARRAY_TASK_ID}]}_trim_log.txt ${INDIR}/${SAMPLE[${SLURM_ARRAY_TASK_ID}]}_R1.fastq.gz ${INDIR}/${SAMPLE[${SLURM_ARRAY_TASK_ID}]}_R2.fastq.gz ${trimDIR}/${SAMPLE[${SLURM_ARRAY_TASK_ID}]}_1.paired_trim.fastq.gz ${trimDIR}/${SAMPLE[${SLURM_ARRAY_TASK_ID}]}_1.unpaired.fastq.gz ${trimDIR}/${SAMPLE[${SLURM_ARRAY_TASK_ID}]}_2.paired_trim.fastq.gz ${trimDIR}/${SAMPLE[${SLURM_ARRAY_TASK_ID}]}_2.unpaired.fastq.gz ILLUMINACLIP:/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10:8:TRUE LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

## removing host genome

bbsplit.sh build=1 threads=8 -Xmx=32g minid=0.95 maxindel=3 bandwidthratio=0.16 bandwidth=12 quickmatch fast minhits=2 pigz unpigz path=/References/bbmap_masked_ref/hg19-OryCun2 in=${trimDIR}/${SAMPLE[${SLURM_ARRAY_TASK_ID}]}_1.paired_trim.fastq.gz in2=${trimDIR}/${SAMPLE[${SLURM_ARRAY_TASK_ID}]}_2.paired_trim.fastq.gz basename=${cleanDIR}/${SAMPLE[${SLURM_ARRAY_TASK_ID}]}.dirty.%.fq.gz outu=${cleanDIR}/${SAMPLE[${SLURM_ARRAY_TASK_ID}]}.clean.fastq.gz

## filtering out broken reads

repair.sh tossbrokenreads  in=${cleanDIR}/${SAMPLE[${SLURM_ARRAY_TASK_ID}]}.clean.fastq.gz out1=${cleanDIR}/${SAMPLE[${SLURM_ARRAY_TASK_ID}]}_1.clean.fastq.gz out2=${cleanDIR}/${SAMPLE[${SLURM_ARRAY_TASK_ID}]}_2.clean.fastq.gz

## alignment with BWA

bwa mem -t 6 ${REF} ${cleanDIR}/${SAMPLE[${SLURM_ARRAY_TASK_ID}]}_1.clean.fastq.gz ${cleanDIR}/${SAMPLE[${SLURM_ARRAY_TASK_ID}]}_2.clean.fastq.gz -R "@RG\\tID:bwa\\tPL:illumina\\tLB:${SAMPLE[${SLURM_ARRAY_TASK_ID}]}_lib\\tSM:${SAMPLE[${SLURM_ARRAY_TASK_ID}]}" | samtools view -Sb -@ 2 -o ${alignDIR}/${SAMPLE[${SLURM_ARRAY_TASK_ID}]}.nic.raw.bam

## sorting and indexing BAM files

samtools sort -@8 ${alignDIR}/${SAMPLE[${SLURM_ARRAY_TASK_ID}]}.nic.raw.bam -o ${alignDIR}/${SAMPLE[${SLURM_ARRAY_TASK_ID}]}.nic.sorted.bam
samtools index -@ 8 ${alignDIR}/${SAMPLE[${SLURM_ARRAY_TASK_ID}]}.nic.sorted.bam

## filtering reads

samtools flagstat -@ 8 ${alignDIR}/${SAMPLE[${SLURM_ARRAY_TASK_ID}]}.nic.sorted.bam > ${alignDIR}/${SAMPLE[${SLURM_ARRAY_TASK_ID}]}.nic.flagstat

samtools view -@ 3 -h -F 268 -f 2 -b ${alignDIR}/${SAMPLE[${SLURM_ARRAY_TASK_ID}]}.nic.sorted.bam | samtools sort -n -@ 3 - | samtools fixmate -O bam - - | samtools sort -@ 1 - > ${alignDIR}/${SAMPLE[${SLURM_ARRAY_TASK_ID}]}.nic.filt.bam

samtools index -@ 6 ${alignDIR}/${SAMPLE[${SLURM_ARRAY_TASK_ID}]}.nic.filt.bam | samtools flagstat ${alignDIR}/${SAMPLE[${SLURM_ARRAY_TASK_ID}]}.nic.filt.bam > ${alignDIR}/${SAMPLE[${SLURM_ARRAY_TASK_ID}]}.nic.filt.flagstat

## marking duplicates

picard  MarkDuplicates \
      I=${alignDIR}/${SAMPLE[${SLURM_ARRAY_TASK_ID}]}.nic.filt.bam O=${alignDIR}/${SAMPLE[${SLURM_ARRAY_TASK_ID}]}.nic.dedup.bam \
      METRICS_FILE=${alignDIR}/${SAMPLE[${SLURM_ARRAY_TASK_ID}]}.nic.dedup.metrics \
      TMP_DIR=${alignTMP} REMOVE_DUPLICATES=true \
      MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=100

## verifying mate-pair information between reads

picard  FixMateInformation INPUT=${alignDIR}/${SAMPLE[${SLURM_ARRAY_TASK_ID}]}.nic.dedup.bam OUTPUT=${alignDIR}/${SAMPLE[${SLURM_ARRAY_TASK_ID}]}.nic.matefixed.bam TMP_DIR=${alignTMP}

picard  BuildBamIndex INPUT=${alignDIR}/${SAMPLE[${SLURM_ARRAY_TASK_ID}]}.nic.matefixed.bam OUTPUT=${alignDIR}/${SAMPLE[${SLURM_ARRAY_TASK_ID}]}.nic.matefixed.bam.bai TMP_DIR=${alignTMP}

## realign indel

gatk -T RealignerTargetCreator \
     -R ${REF}  -I ${alignDIR}/${SAMPLE[${SLURM_ARRAY_TASK_ID}]}.nic.matefixed.bam \
     -o ${alignDIR}/${SAMPLE[${SLURM_ARRAY_TASK_ID}]}.nic.realigner.intervals

gatk -T IndelRealigner \
     -R ${REF} -I ${alignDIR}/${SAMPLE[${SLURM_ARRAY_TASK_ID}]}.nic.matefixed.bam \
     -targetIntervals ${alignDIR}/${SAMPLE[${SLURM_ARRAY_TASK_ID}]}.nic.realigner.intervals \
     -o ${alignDIR}/${SAMPLE[${SLURM_ARRAY_TASK_ID}]}.nic.realn.bam

## removing excessive mismatches

bamutils filter ${alignDIR}/${SAMPLE[${SLURM_ARRAY_TASK_ID}]}.nic.realn.bam ${filterDIR}/${SAMPLE[${SLURM_ARRAY_TASK_ID}]}.nic.mm.bam -failed ${filterDIR}/${SAMPLE[${SLURM_ARRAY_TASK_ID}]}.nic.mm.txt -mismatch 5

## filtering short reads

bamutils filter ${filterDIR}/${SAMPLE[${SLURM_ARRAY_TASK_ID}]}.nic.mm.bam ${filterDIR}/${SAMPLE[${SLURM_ARRAY_TASK_ID}]}.nic.minlen.bam -minlen 35

## removing soft-clips

bamutils removeclipping ${filterDIR}/${SAMPLE[${SLURM_ARRAY_TASK_ID}]}.nic.minlen.bam ${filterDIR}/${SAMPLE[${SLURM_ARRAY_TASK_ID}]}.nic.scf.tmp.bam

samtools view -@ 4 ${filterDIR}/${SAMPLE[${SLURM_ARRAY_TASK_ID}]}.nic.scf.tmp.bam | grep 'ZC:f:' | awk '{for (i=1;i<=NF;i++){if ($sample ~/ZC:f:/) {print $1, $sample}}}' | sed 's/ZC:f://' | awk -v MAX_SOFTCLIP="0.05" '$2 > MAX_SOFTCLIP {print $1}' > ${filterDIR}/${SAMPLE[${SLURM_ARRAY_TASK_ID}]}.read_names_to_remove_highSoftClip.txt.tmp

cat ${filterDIR}/${SAMPLE[${SLURM_ARRAY_TASK_ID}]}.read_names_to_remove_highSoftClip.txt.tmp | cut -f 1 | sort -T ${alignTMP} --parallel=4 | uniq > ${filterDIR}/${SAMPLE[${SLURM_ARRAY_TASK_ID}]}.read_names_to_remove_highSoftClip.txt

picard FilterSamReads I=${filterDIR}/${SAMPLE[${SLURM_ARRAY_TASK_ID}]}.nic.minlen.bam O=${filterDIR}/${SAMPLE[${SLURM_ARRAY_TASK_ID}]}.nic.scf.filtOut.bam READ_LIST_FILE=${filterDIR}/${SAMPLE[${SLURM_ARRAY_TASK_ID}]}.read_names_to_remove_highSoftClip.txt FILTER=includeReadList

picard FilterSamReads I=${filterDIR}/${SAMPLE[${SLURM_ARRAY_TASK_ID}]}.nic.minlen.bam O=${filterDIR}/${SAMPLE[${SLURM_ARRAY_TASK_ID}]}.nic.scf.bam READ_LIST_FILE=${filterDIR}/${SAMPLE[${SLURM_ARRAY_TASK_ID}]}.read_names_to_remove_highSoftClip.txt FILTER=excludeReadList

## removing chimeric alignment

samtools view -@ 4 -h -f 2048 -b ${filterDIR}/${SAMPLE[${SLURM_ARRAY_TASK_ID}]}.nic.scf.bam > ${filterDIR}/${SAMPLE[${SLURM_ARRAY_TASK_ID}]}.nic.sup.filtOut.bam & samtools view -@ 4 -h -F 2048 -b ${filterDIR}/${SAMPLE[${SLURM_ARRAY_TASK_ID}]}.nic.scf.bam > ${filterDIR}/${SAMPLE[${SLURM_ARRAY_TASK_ID}]}.nic.sup.bam

## removing hard-clips

samtools view -@ 4 -h ${filterDIR}/${SAMPLE[${SLURM_ARRAY_TASK_ID}]}.nic.sup.bam | awk '$6 ~/H/{print}' | samtools view -Sbh -  > ${filterDIR}/${SAMPLE[${SLURM_ARRAY_TASK_ID}]}.nic.hcf.filtOut.bam & samtools view -@ 4 -h ${filterDIR}/${SAMPLE[${SLURM_ARRAY_TASK_ID}]}.nic.sup.bam | awk ' 6 !~ /H/{print} ' | samtools view -bh > ${filterDIR}/${SAMPLE[${SLURM_ARRAY_TASK_ID}]}.nic.hcf.bam

## MAPQ10 filter

samtools view -@ 4 -h ${filterDIR}/${SAMPLE[${SLURM_ARRAY_TASK_ID}]}.nic.hcf.bam | awk -v var="10" '$5 < var || $1 ~/^@/' | samtools view -b - > ${filterDIR}/${SAMPLE[${SLURM_ARRAY_TASK_ID}]}.nic.MAPQ10.filtOut.bam & samtools view -@ 4 -h -bq 10 ${filterDIR}/${SAMPLE[${SLURM_ARRAY_TASK_ID}]}.nic.hcf.bam > ${filterDIR}/${SAMPLE[${SLURM_ARRAY_TASK_ID}]}.nic.MAPQ10.bam

## MAPQ40 filter

samtools view -@ 4 -h ${filterDIR}/${SAMPLE[${SLURM_ARRAY_TASK_ID}]}.nic.MAPQ10.bam | awk -v var="40" '$5 < var || $1 ~/^@/' | samtools view -b - > ${filterDIR}/${SAMPLE[${SLURM_ARRAY_TASK_ID}]}.nic.MAPQ40.filtOut.bam & samtools view -@ 4 -h -bq 40 ${filterDIR}/${SAMPLE[${SLURM_ARRAY_TASK_ID}]}.nic.MAPQ10.bam > ${filterDIR}/${SAMPLE[${SLURM_ARRAY_TASK_ID}]}.nic.MAPQ40.bam

## removing singletone

samtools view -@ 4 ${filterDIR}/${SAMPLE[${SLURM_ARRAY_TASK_ID}]}.nic.MAPQ40.bam | cut -f 1 | sort -T ${alignTMP}/ --parallel=4 | uniq -u > ${filterDIR}/${SAMPLE[${SLURM_ARRAY_TASK_ID}]}.singleAfterFilt.txt

picard FilterSamReads I=${filterDIR}/${SAMPLE[${SLURM_ARRAY_TASK_ID}]}.nic.MAPQ40.bam O=${filterDIR}/${SAMPLE[${SLURM_ARRAY_TASK_ID}]}.nic.singletons.filtOut.bam READ_LIST_FILE=${filterDIR}/${SAMPLE[${SLURM_ARRAY_TASK_ID}]}.singleAfterFilt.txt FILTER=includeReadList

picard FilterSamReads I=${filterDIR}/${SAMPLE[${SLURM_ARRAY_TASK_ID}]}.nic.MAPQ40.bam O=${filterDIR}/${SAMPLE[${SLURM_ARRAY_TASK_ID}]}.nic.MAPQ40.singletons.final.bam READ_LIST_FILE=${filterDIR}/${SAMPLE[${SLURM_ARRAY_TASK_ID}]}.singleAfterFilt.txt FILTER=excludeReadList

## variant calling

gatk -T HaplotypeCaller \
     -R ${REF} \
     -I ${bamDIR}/${SAMPLE[${SLURM_ARRAY_TASK_ID}]}.nic.MAPQ40.singletons.final.bam \
     -ploidy 1 -stand_call_conf 30 --emitRefConfidence GVCF \
     -o ${gvcfDIR}/${SAMPLE[${SLURM_ARRAY_TASK_ID}]}_HaploCaller.raw.g.vcf.gz \
     --unsafe

gatk -T CombineGVCFs -R ${REF} --variant ${DIR}/TPA_gvcfList.txt -o ${DIR}/TPA_HaploCaller_joint.raw.g.vcf

gatk -T GenotypeGVCFs -R ${REF} --variant ${DIR}/TPA_HaploCaller_joint.raw.g.vcf -o ${DIR}/TPA_HaploCaller_joint.raw.vcf

gatk -T VariantFiltration -R ${REF} -V ${DIR}/TPA_HaploCaller_joint.raw.vcf \
     --filterExpression "QD < 2.0" --filterName "QD" \
     --filterExpression "MQ < 40.0" --filterName "MQ" \
     --filterExpression "FS > 60.0" --filterName "FS" \
     --filterExpression "SOR > 3.0" --filterName "SORhigh" \
     --filterExpression "MQRankSum <  -12.5" --filterName "MQRankSum" \
     --filterExpression "ReadPosRankSum < -8.0 " --filterName "ReadPosRankSum" \ 
     --filterExpression "DP < 3.0" --filterName "LowCoverage" \
     --logging_level ERROR \
     -o ${DIR}/TPA_HaploCaller_joint.qual.vcf

gatk -T SelectVariants -R ${REF} \
     -V ${DIR}/TPA_HaploCaller_joint.qual.vcf \
     -o ${DIR}/TPA_HaploCaller_joint.pass.vcf \
     -select "vc.isNotFiltered()"

