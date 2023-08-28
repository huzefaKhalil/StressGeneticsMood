#!/bin/bash

# set the directories
inDir="../qc/"
inFile="merge_6"

outDir="../imputationFiles/"
outFile="mergedForImputation"

# let's create the directory
mkdir -p $outDir

# now, set the files
inFile=$inDir$inFile
outFile=$outDir$outFile

echo "Starting plink in directory $inFile"

perlscript="HRC-1000G-check-bim.pl"
hrcref="../../ref/HRC.r1-1.GRCh37.wgs.mac5.sites.tab"
#onekgref="../ref/1000GP_Phase3_combined.legend"

#plink --bfile ${inFile} --freq --out ${outFile}_frq

# you could replace $onekgref to $hrcref if using HRC as the imputation panel
#perl $perlscript -g -b ${file_name}.bim -f ${file_name}_frq.frq -r $onekgref -v 
#perl $perlscript -h -b ${inFile}.bim -f ${outFile}_frq.frq -r $hrcref -v

# move files to outDir
mv *-HRC.txt $outDir

# Modify the Run plink file generated to read and write to the correct directories
sed -i "s| merge_6 | $inFile |g" Run-plink.sh
sed -i "s|Exclude-merge_6-HRC.txt|${outDir}Exclude-merge_6-HRC.txt|g" Run-plink.sh
sed -i "s|Chromosome-merge_6-HRC.txt|${outDir}Chromosome-merge_6-HRC.txt|g" Run-plink.sh
sed -i "s|Position-merge_6-HRC.txt|${outDir}Position-merge_6-HRC.txt|g" Run-plink.sh
sed -i "s|Strand-Flip-merge_6-HRC.txt|${outDir}Strand-Flip-merge_6-HRC.txt|g" Run-plink.sh
sed -i "s|Force-Allele1-merge_6-HRC.txt|${outDir}Force-Allele1-merge_6-HRC.txt|g" Run-plink.sh
sed -i "s|merge_6-updated|${outDir}merge_6|g" Run-plink.sh

# Run plink and get vcfs for all chrs
sh Run-plink.sh


#sed -i -e 's/plink/plink-1.9/g' Run-plink.sh
#sed -i '/23/d' Run-plink.sh

