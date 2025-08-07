#!/bin/bash

# TO run this script:
# sbatch  -p gpu --mem=24g  --gres=gpu:k80:1 hiccups.sh

module load juicer
module load CUDA/10

HIC_MATRIX=$1
OUTPUT_DIR=$2

juicer_tools  hiccups -m 512 -r 5000,10000,25000 -k KR \
    -f .1,.1,.1 -p 4,2,1 -i 7,5,3 -t 0.02,1.5,1.75,2 -d 20000,20000,50000 \
    $HIC_MATRIX $OUTPUT_DIR

# use juicer_tools48g if more memory is required
