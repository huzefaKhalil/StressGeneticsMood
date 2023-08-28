# StressGeneticsMood
This repository contains the code for the analysis presented in the paper "Stress, Genetics and Mood: Impact of COVID-19 on a College Freshman Sample"

## Genetics
The bash scripts for manipulating the genetic data require the following software and resources:
* [plink][https://www.cog-genomics.org/plink2/]
* [QC R scripts from Marees et al (2019)][https://github.com/MareesAT/GWA_tutorial/]
* [Check VCF][https://github.com/zhanxw/checkVCF]
* [HRC or 1000G Imputation preparation and checking][https://www.well.ox.ac.uk/~wrayner/tools/]
* [VCF Tools][https://vcftools.github.io/perl_examples.html]
* [Michigan Imputation Server][https://imputationserver.sph.umich.edu/index.html]
* [BCF tools and HTS lib][http://www.htslib.org/download/]
* [PRSice2][https://choishingwan.github.io/PRSice/]

The QC folder contains the scripts to perform quality control on the genetic data and split the data by chromosome in order to get it ready for imputation.
The PostImputation folder contains scripts to combine data from various cohorts after imputation.
The PRSice folder contains the script used to compute the PRS using the [PRSice2][https://choishingwan.github.io/PRSice/] software.
Creating the Polygenic Risk Score requires the GWAS summary statistics from [Howard et al (2019)][https://doi.org/10.1038/s41593-018-0326-7]. The summary statistics for all variants was used.

## R Code for the figures
The R code used to generate all the figures in the paper is provided here. The code needs the survey data which is available upon request. The R utility functions used to access the survey data will also be made available upon request.
