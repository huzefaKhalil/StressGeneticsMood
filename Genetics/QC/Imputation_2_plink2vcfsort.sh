#!/bin/bash

# set the directories
inDir="../imputationFiles/"
inFile="merge_6"

outDir="../imputationFiles/"
outFile="merge_6"

# let's create the directory
mkdir -p $outDir

# now, set the files
inFile=$inDir$inFile
outFile=$outDir$outFile

####### convert to vcf and sort
for x in {1..22}
do
	echo ${x}
 	vcf-sort ${inFile}-chr${x}.vcf | bgzip -c > ${outFile}-chr${x}.vcf.gz
done
