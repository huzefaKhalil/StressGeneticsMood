#!/bin/bash

gwasFile="../ref/imputedGwasRefFileV4.txt"
targetFile="../combinedDataP4/imputedOutput2/"

outDir="../combinedDataP4/clumpedPRS/"
outFile="${outDir}prsRes"

phenotypeFile="../combinedDataP4/combined.phenotype"

mkdir -p $outDir


# first we run it to make the *.valid file... then we run it to make the PRS
for i in {1..22}
do

/usr/bin/Rscript PRSice.R --dir PRSiceV2 \
        --prsice PRSice_linux \
        --base $gwasFile \
        --target ${targetFile}chr${i}bgen,${targetFile}chr${i}bgen.sample \
        --snp SNP \
        --A1 A1 \
        --A2 A2 \
        --bp BP \
        --chr CHR \
        --stat LogOR\
        --beta \
        --pvalue P \
        --base-info Rsq:0.9 \
        --base-maf MAF:0.1 \
        --type bgen \
        --pheno-file $phenotypeFile \
        --ignore-fid \
        --pheno-col Phen1 \
        --thread max \
        --score sum \
        --binary-target F \
        --fastscore \
        --model add \
        --all-score \
        --out  ${outDir}chr${i} \
        --bar-levels 0.0001,0.00001 \
        --lower 0.00001 \
        --interval 0.00001 \
        --upper 1 \
        --missing SET_ZERO \
        --extract ../combinedDataP4/clumpedPRS/snp.valid
done

## now run it as it will actually work... there must be a better way to do this
#for i in {1..22}
#do
#
#Rscript PRSice.R --dir PRSiceV2 \
#	--prsice PRSice_linux \
#	--base $gwasFile \
#	--target ${targetFile}chr${i}bgen,${targetFile}chr${i}bgen.sample \
#	--snp SNP \
#	--A1 A1 \
#	--A2 A2 \
#	--bp pos \
#	--chr seqnames \
#	--stat LogOR\
#	--beta \
#	--pvalue P \
#	--base-info Rsq,0.9 \
#	--base-maf MAF,0.1 \
#	--type bgen \
#	--pheno-file $phenotypeFile \
#	--ignore-fid \
#	--pheno-col Phen1 \
#	--thread 1 \
#	--score sum \
#	--binary-target F \
#	--fastscore \
#	--model add \
#	--all-score \
#	--no-clump \
#	--out  ${outDir}chr${i} \
#	--bar-levels 0.0001,0.00001 \
#	--lower 0.00001 \
#	--interval 0.00001 \
#	--upper 1 \
#	--extract ${outDir}chr${i}.valid \
#	--missing SET_ZERO
#
#done
