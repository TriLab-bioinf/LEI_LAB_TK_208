#!/usr/bin/bash

# check_potential_restriction_enzyme.sh <gzipped R1 fastq file> <number of bases to check from 5'> <number of reads to evaluate>

FASTQ=$1
LEN=$2
NUM_SEQ=$3

for seq in  $(gunzip -c ${FASTQ} | \
    head -${NUM_SEQ} | \
    paste - - - - | \
    cut -f 2 ) 
do 
    echo ${seq:0:${LEN} }; 
done | \
    sort |\
    uniq -c | \
    sort -k1n
