#!/bin/bash

# This file combines the files post imputation.
# It should be run after unzipping the files

# first the directories
data1="../data1/imputed1/"
data2="../data2/imputed1/"
data3="../data3/imputed1/"
data4="../data4/imputed1/"

outData="../combinedDataP4/imputedOutput/"

mkdir -p $outData

# 1. Merge using bcftools

for i in {1..22}
do
    echo "Starting chr $i"

    # first create the indices

    bcftools index -f "${data1}chr${i}.dose.vcf.gz"
    bcftools index -f "${data2}chr${i}.dose.vcf.gz"
    bcftools index -f "${data3}chr${i}.dose.vcf.gz"
    bcftools index -f "${data4}chr${i}.dose.vcf.gz"

    # next, merge the files
    bcftools merge -i AF:avg,MAF:avg,R2:avg -m id "${data1}chr${i}.dose.vcf.gz" "${data2}chr${i}.dose.vcf.gz" "${data3}chr${i}.dose.vcf.gz" "${data4}chr${i}.dose.vcf.gz" -Oz -o "${outData}chr${i}.merged.vcf.gz"

done
