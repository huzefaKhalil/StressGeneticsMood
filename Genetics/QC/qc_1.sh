#!/bin/bash 

inDir="../input/"
outDir="../qc/"

mkdir -p $outDir

inFile="merge"

outFile=${outDir}${inFile}
inFile=${inDir}${inFile}

####################################################
# STEP 0: add sex to the ped file - Don't need to do this here as its done in the merged file. 
#awk 'FNR==NR{a[NR]=$5;next}{$5=a[FNR]}1' ${inFile}.fam ${inFile}.ped > tmp.file
#mv tmp.file ${inFile}.ped

####################################################
# STEP 1: check missingness and generate plots
plink --bfile ${inFile} --missing --out $outFile

/usr/bin/Rscript --no-save hist_miss.R ${outFile}.imiss ${outFile}.lmiss ${outDir}

####################################################
# STEP 2: remove individuals with high missingness.. They recommend using a less stringent filter and then following it with a more stringent filter. Even the more stringent filter is more relaxed than Srijan's
# Also, they filter on SNP missingness and indivudal missingness. Srijan first filters on individual missingness and then on SNPs much later

# filter SNPs at 0.2
plink --bfile ${inFile} --geno 0.1 --make-bed --out ${outFile}_1

# filter individuals at 0.2
plink --bfile ${outFile}_1 --mind 0.05 --make-bed --out ${outFile}_2

####################################################
# STEP 3: sexcheck
plink --bfile ${outFile}_2 --check-sex --out ${outFile}

# visualize the sex check
/usr/bin/Rscript --no-save sex_check.R ${outFile}.sexcheck ${outDir}

# remove individuals with problematic sex
grep "PROBLEM" ${outFile}.sexcheck | awk '{print$1,$2}' > ${outDir}sex_discrepency.txt
plink -bfile ${outFile}_2 --remove ${outDir}sex_discrepency.txt --make-bed --out ${outFile}_3

# impute sex... (not sure why this is necessary). NOT RUNNING IT
#plink -bfile ${outFile}_5 --impiute-sex --make-bed --out ${outFile}_6

####################################################
# STEP 4: select autosomal SNPs only an filter out SNPs with low minor allele frequency (MAF)

# select autosomal SNPs only, ie from chr 1 to 22
awk '{ if ($1 >= 1 && $1 <= 22) print $2 }' ${outFile}_3.bim > ${outDir}snp_1_22.txt
plink --bfile ${outFile}_3 --extract ${outDir}snp_1_22.txt --make-bed --out ${outFile}_4

# generat plot of MAF distribution
plink --bfile ${outFile}_4 --freq --out ${outDir}MAF_check

#visualize it
/usr/bin/Rscript --no-save MAF_check.R ${outDir}MAF_check. ${outDir}

# remove SNPs with low MAF... major point of diversion. Srijan's MAF filtering crieria is VERY small. 0.005 vs what they recommend here of 0.05. I'll go midway with 0.01.
plink --bfile ${outFile}_4 --maf 0.005 --make-bed --out ${outFile}_5

####################################################
# STEP 5: Delete SNPs not in the Hardy-WEinberg equilibrium (HWE)

plink --bfile ${outFile}_5 --hardy --out ${outFile}

# select SNPs with HWE p-value below 0.00001
awk '{ if ($9 < 0.0001) print $0 }' ${outFile}.hwe > ${outFile}zoomhwe.hwe
/usr/bin/Rscript --no-save hwe.R ${outFile}.hwe ${outFile}zoomhwe.hwe ${outDir}

# now delete them. We don'd have cae / controls, so filter all at 1e-10
# this is again a departure from Yu's pipeline as they filter at 1e-6
#plink --bfile ${outFile}_7 --hwe 1e-6 --make-bed --out ${outFile}_8
plink --bfile ${outFile}_5 --hwe 1e-10 --hwe-all --make-bed --out ${outFile}_6

############################################################
# STEP 6: Heterozygosity and LD Pruning

# so first with heterzygosity rate. Remove individuals with a heterozygosity rate of more than 3 sd from the mean
# heterozygosity: carrying two different alleles of a SNP. Low indicates inbreeding, high indicates low samples quality
# exclude inversion regions (high LD regions).
# Linkage Disequilibrium (LD) a measure of non-random association between alleles at different loci in the same chromosome.
# SNPs are in LD when the frequency of association of their alleles is higher than expected under random assortment.
# To calculate LD, we use the following parameters: window size=50, shift size=5, pairwise SNP-SNP r^2 metric=0.2
# Yu's parameters are 100, 25, 0.5. Given the small sample sizes we are dealing with, Yu's parameters, particularly for r^2
# make more sense. More SNPs would be excluded with r^2=0.2

plink --bfile ${outFile}_6 --indep-pairwise 50 5 0.5 --out ${outDir}indepSNP

# get the pruned data set
plink --bfile ${outFile}_6 --extract ${outDir}indepSNP.prune.in --het --out ${outDir}R_hetCheck

# plot the heterogygosity rate
RScript --no-save heterogygosity_rate.R ${outDir}R_hetCheck ${outDir}

# list the individuals more than 3 sd away from the heterozygosity rate mean.
/usr/bin/Rscript --no-save heterozygosity_outliers.R ${outDir}R_hetCheck ${outDir}

# need to exclude these individuals from the analysis.
sed 's/"// g' ${outDir}fail-het-qc.txt | awk '{print $1, $2}' > ${outDir}het-fail-ind.txt

plink --bfile ${outFile}_6 --remove ${outDir}het-fail-ind.txt --make-bed --out ${outFile}_7

############################################################
# STEP 7: Relatedness

# This step isn't there in Yu's workflow, because it shouldn't really be necessary.

plink --bfile ${outFile}_7 --extract ${outDir}indepSNP.prune.in --genome --min 0.2 --out ${outDir}pihat_min0.2

# visualize relations. But there will be none! or should be none!
awk '{ if ($8 > 0.9) print $0}' ${outDir}pihat_min0.2.genome > ${outDir}zoom_pihat.genome

/usr/bin/Rscript --no-save relatedness.R ${outDir}pihat_min0.2.genome ${outDir}zoom_pihat.genome ${outDir}

# we shouldn't have any parent offspring relationships! but let's see if anything pops up
awk '{ if ($8 > 0.9) print $0 }' ${outDir}pihat_min0.2.genome > ${outDir}zoom_pihat.genome

# visualize
/usr/bin/Rscript --no-save relatedness.R ${outDir}pihat_min0.2.genome ${outDir}zoom_pihat.genome ${outDir}


#if there is anyone with a piHat more than 0.2, remove them!


