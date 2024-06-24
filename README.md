# Bayesian Bootstrap Analyses

This repository contains the code of the "Bayesian bootstrap" analyses.

It is organized in 6 subdirectories corresponding to the 6 use-cases:

- 01_simple_example: simple example of high supports of null branches and low support of short branches.
- 02_covid: Analysis of SARS-CoV-2 dataset.
- 03_ebola: Analysis of the EBOV dataset.
- 04_RVFV: Analysis of the Rift Valley Fever Virus dataset.
- 05_real_datasets: Analysis of 300+100 datasets from https://cme.h-its.org/exelixis/material/raxml_adaptive_data.tar.gz .
- 06_simulations: Aalysis of simulated datasets with different levels of homoplasy.

To run the different analyses, a run.sh script is present in each folder, describing the steps. It generally involves running a nextflow workflow.

The required software are: Java and Apptainer/Singularity.

The results are given in the release of this repository.
