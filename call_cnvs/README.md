# Calling GYP CNVs

## Introduction
We have implemented a Hidden Markov model (HMM) to identify copy number variation from sequence coverage. The approach has been applied to call copy number variants in the 360 kb segmental duplication on chromosome 4 that encodes the glycophorin genes *GYPE*, *GYPB* and *GYPA*. Coverage data was pre-processed to exclude sites and windows with poor mappability (see description of example data). Multiple parameter choices are hard coded and this method has not been applied in other regions. This serves to document the implementation in the paper but it has not been developed for general use. Two R scripts are provided - one to run the HMM and output copy number state paths, and one to plot the results.

## Running the HMM
The HMM R script should be run with command line arguments of the input file containing coverage data, and an optional output file prefix.  If no output file prefix is supplied, "paths" will be used by default. This script requires the "msm" package to be installed.

```
Rscript --vanilla gypcnvs_hmm.R inputfile [outprefix]
```

### Input data format
The input file is a matrix of the mean coverage for each individual (rows) in each bin (columns) from a region on a single chromosome. The column names are the start position of each bin and the row names are individual IDs. NA is used to indicate a bin that has been excluded.

The script reads the input data, runs the HMM, and outputs a matrix of copy number state for each individual (rows) in each bin (columns) to the file *"outprefix.out.txt"*.

## Plotting the results
Similarly, the plotting script should be run with command line arguments of the input file containing copy number state paths, and an optional output file prefix. If no output file prefix is supplied, "paths" will be used by default.

```
Rscript --vanilla plot_cnvs.R inputfile [outprefix]
```

This will read the file of copy number state paths, and produce a plot for individuals with copy number variants in the file *"outprefix.plot.pdf"*.

## Example data
The provided example input file, "inds200.cov.txt", contains coverage data across the glycophorin region for 200 1000 Genomes Phase 3 individuals subset from the study. This includes 98 individuals with a copy number variant and 102 individuals with no copy number variant. Coverage was calculated using bedtools genomecov and then averaged over 1600 bp windows.  Because of the high sequence identity between the genes and intergenic sequences, sites with low mappability (mappability score \<0.9 where mappability of a site is the mean value of the CRG mappability track for all 100-mers overlapping that site) and bins with fewer than 25% such mappable sites were excluded.  To run the scripts on this example data, you would run these two lines:
```
Rscript --vanilla gypcnvs_hmm.R inds200.cov.txt
Rscript --vanilla plot_cnvs.R paths.out.txt
```
The output files should match the two files in the directory "sample_output_files".
