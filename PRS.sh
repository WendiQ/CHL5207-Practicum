#!/bin/bash

cd /hpf/largeprojects/struglis/wendi/PRS

module load plink/2.0
module load plink/1.90b6.21


plink \
    --bfile Merged-TOPMed-7arrays-10XG-PacBioHIFI-N5467-n7217886-postQC \
    --score PRS_weights_hg38_common.txt 1 2 3 header sum \
    --keep subjects_keep.txt \
    --out PRS_gen_cf


