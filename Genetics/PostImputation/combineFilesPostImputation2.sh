#!/bin/bash

# Second file, this makes the info with the ID, and MAF and Rsq

# so, let's make the info file

inDir="../combinedDataP4/imputedOutput/"

outDir="../combinedDataP4/imputedOutput2/"

mkdir -p $outDir

outInfoFile="${outDir}allRsqMaf.info"

echo -e 'SNP\tMAF\tRsq' > $outInfoFile

for i in {1..22}
do
    bcftools query -f '%ID\t%INFO/MAF\t%INFO/R2\n' "${inDir}chr${i}.merged.vcf.gz" >> $outInfoFile
done

# now generate indices
echo "Generating indices"

for i in {1..22}
do
    tabix -p vcf "${inDir}chr${i}.merged.vcf.gz"
done

# now convert to bgen
echo "Converting to bgen"
for i in {1..22}
do
    /garage/akil_lab/plink/plink2 --vcf "${inDir}chr${i}.merged.vcf.gz" dosage=DS --export bgen-1.1 --out "${outDir}chr${i}bgen"
done

# now add sex to the sample file
#sexFile="example.sample"
#echo "Adding sex"
#for i in {1..22}
#do
#    awk 'FNR==NR{a[NR]=$4;next}{$4=a[FNR]}1' $sexFile "${outDir}chr${i}bgen.sample" > tmp.file
#    mv tmp.file "${outDir}chr${i}bgen.sample"
#done

#Rscript PostImputationRSidConversion.R
