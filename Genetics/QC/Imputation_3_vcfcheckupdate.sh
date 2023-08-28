#!/bin/bash

# set the directories
inDir="../imputationFiles/"
inFile="merge_6"

outDir="../imputationFilesForServer/"
outFile="merge_6"

refDir="../../ref/"
refFile=${refDir}human_g1k_v37.fasta

# let's create the directory
mkdir -p $outDir

# now, set the files
inFile=$inDir$inFile
outFile=$outDir$outFile

for x in {1..22}
do
 	echo ${x}
 	python2.7 checkVCF.py -r ${refFile} -o ${outFile}-chr${x} ${inFile}-chr${x}.vcf.gz
done

cat ${outDir}*.log > ${outDir}all.logfile
