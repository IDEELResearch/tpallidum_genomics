#!/bin/bash
#SBATCH -n 1
#SBATCH -t 6:00:00
#SBATCH --mem 64Gb

projdir=/pipeline/postfilter

$projdir/1.remove_too_much_mismatches.sh
#$projdir/2.remove_short_alignment.sh
#$projdir/3.too_much_softclipped.sh
#$projdir/4.sulementary_chimeric_alignment.sh
#$projdir/5.hardclipped.sh
#$projdir/6.MAPQ10.sh
#$projdir/7.MAPQ40.sh
#$projdir/8.singletons.sh

#

#snakemake -s do-alignment-ngmlr.py --cluster "sbatch -n5 -t 4-00:00:00 --mem 7Gb " -j 8
