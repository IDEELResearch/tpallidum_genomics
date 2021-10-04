#!/usr/bin/bash

#extract multiple sequences from a large fasta file
module load seqkit

while read line; do less batch2_ss14_ompeome_all.fa | seqkit grep -p $line; done <TP0126-OmpW1.txt >>TP0126-OmpW1.fasta
while read line; do less batch2_ss14_ompeome_all.fa | seqkit grep -p $line; done <TP0326-BamA.txt >>TP0326-BamA.fasta
while read line; do less batch2_ss14_ompeome_all.fa | seqkit grep -p $line; done <TP0479-OprG.txt >>TP0479-OprG.fasta
while read line; do less batch2_ss14_ompeome_all.fa | seqkit grep -p $line; done <TP0515-LptD.txt >>TP0515-LptD.fasta
while read line; do less batch2_ss14_ompeome_all.fa | seqkit grep -p $line; done <TP0548-FadL.txt >>TP0548-FadL.fasta
while read line; do less batch2_ss14_ompeome_all.fa | seqkit grep -p $line; done <TP0698-New_OMP.txt >>TP0698-New_OMP.fasta
while read line; do less batch2_ss14_ompeome_all.fa | seqkit grep -p $line; done <TP0705-mrcA_gene.txt >>TP0705-mrcA_gene.fasta
while read line; do less batch2_ss14_ompeome_all.fa | seqkit grep -p $line; done <TP0733-OmpW2.txt >>TP0733-OmpW2.fasta
while read line; do less batch2_ss14_ompeome_all.fa | seqkit grep -p $line; done <TP0751-laminin-binding_protein.txt >>TP0751-laminin-binding_protein.fasta
while read line; do less batch2_ss14_ompeome_all.fa | seqkit grep -p $line; done <TP0856-FadL.txt >>TP0856-FadL.fasta
while read line; do less batch2_ss14_ompeome_all.fa | seqkit grep -p $line; done <TP0858-FadL.txt >>TP0858-FadL.fasta
while read line; do less batch2_ss14_ompeome_all.fa | seqkit grep -p $line; done <TP0859-FadL.txt >>TP0859-FadL.fasta
while read line; do less batch2_ss14_ompeome_all.fa | seqkit grep -p $line; done <TP0865-FadL.txt >>TP0865-FadL.fasta
while read line; do less batch2_ss14_ompeome_all.fa | seqkit grep -p $line; done <TP0966-OprJ.txt >>TP0966-OprJ.fasta
while read line; do less batch2_ss14_ompeome_all.fa | seqkit grep -p $line; done <TP0967-OprN.txt >>TP0967-OprN.fasta
while read line; do less batch2_ss14_ompeome_all.fa | seqkit grep -p $line; done <TP0968-New_OMP.txt >>TP0968-New_OMP.fasta
while read line; do less batch2_ss14_ompeome_all.fa | seqkit grep -p $line; done <TP0969-TolC.txt >>TP0969-TolC.fasta

