#!/bin/bash

#SBATCH -n 4
#SBATCH -t 2:00:00
#SBATCH --mem 24Gb
#SBATCH -o %A_pca.out


####################
module load vcftools/0.1.15
module load plink/1.90b6.21
Module load eigensoft/6.1.4
####################

# Keep SNPs

vcftools --gzvcf {FILE_NAME}.final.vcf.gz --out {FILE_NAME}.final.snp --recode --remove-indels --stdout | gzip -c > {FILE_NAME}.final.snp.vcf.gz

# Keep balletic

vcftools --gzvcf {FILE_NAME}.final.snp.vcf.gz --out {FILE_NAME}.final.snp.biallelic --recode --max-alleles 2 --min-alleles 2 --stdout | gzip -c > {FILE_NAME}.final.snp.biallelic.vcf.gz

# Keep SNPs with minor allele frequency (MAF) greater than 0.05

vcftools --gzvcf {FILE_NAME}.final.snp.biallelic.vcf.gz --out {FILE_NAME}.final.snp.biallelic.maf0.05 --recode --maf 0.05 --stdout | gzip -c > {FILE_NAME}.final.snp.biallelic.maf0.05.vcf.gz

# Filter out SNPs with excessive missingness

vcftools --gzvcf {FILE_NAME}.final.snp.biallelic.maf0.05.vcf.gz --out {FILE_NAME}.final.snp.biallelic.maf0.05.miss0.8 --max-missing 0.8 --recode --stdout | gzip -c > {FILE_NAME}.final.snp.biallelic.maf0.05.miss0.8.vcf.gz

# Converting VCF to PLINK file format cell 

plink --vcf {FILE_NAME}.final.snp.biallelic.maf0.05.miss0.8.vcf.gz --make-bed --keep-allele-order --out {FILE_NAME}.final.snp.biallelic.maf0.05.miss0.8 --allow-extra-chr --double-id --vcf-half-call haploid --vcf-idspace-to "_"

# Convert chromosome code to plink format

sed 's/CP004010.2/1/' {FILE_NAME}.final.snp.biallelic.maf0.05.miss0.8.bim > {FILE_NAME}.final.snp.biallelic.maf0.05.miss0.8.tmp
mv {FILE_NAME}.final.snp.biallelic.maf0.05.miss0.8.tmp {FILE_NAME}.final.snp.biallelic.maf0.05.miss0.8.bim

# Convert PLINK BED to PED format

plink --bfile {FILE_NAME}.final.snp.biallelic.maf0.05.miss0.8 --recode --keep-allele-order {FILE_NAME}.final.snp.biallelic.maf0.05.miss0.8

# Convert PLINK to EIGENSOFT format

convertf -p par.u19.PED.EIGENSTRAT 

# Running smartpca

smartpca -p par.u19.smartpca

