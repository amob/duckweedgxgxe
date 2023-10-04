# duckweedgxgxe

This repository contains the input data and analysis scripts for manuscript
"Evolutionary consequences of microbiomes for hosts:impacts on host fitness, traits, and heritability"

Microbiome data was first processed in QIIME2
The input datafiles are in folder "qiimeinputsGxGxE"
The scripts to process these files are in folder "qiimescriptsGxGxE"
*Note that file sepp-refs-gg-13-8.qza is not included. This file was downloaded from https://docs.qiime2.org/2021.4/data-resources/  november 2020.
*Note also that directories must be changed to fit the users system, if re-running
Outputs of these scripts that are used in R scripts are in folder "R inputs"

Experimental data and outputs of QIIME2 were processed with R scripts
The inputs to these scripts are located in "R inputs" which also includes some intermediate files for convenience (see script comments)
The main analysis script is "duckweed GxGxE analysis.R"
*Note that "R inputs" also contains instructions for inspecting the provenance for QIIME2 files.
*Note also that some directory paths may need to be changed if re-running scripts.
