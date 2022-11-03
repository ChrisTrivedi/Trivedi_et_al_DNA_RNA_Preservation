# Trivedi, et al. (2022) DNA/RNA Preservation - code repository

This is the repository for the manuscript ["DNA/RNA Preservation in Glacial Snow and Ice Samples"](https://www.frontiersin.org/articles/10.3389/fmicb.2022.894893/full). The raw sequencing data can be obtained from the Sequence Read Archive at [NCBI](https://www.ncbi.nlm.nih.gov/) under BioProject PRJNA657180 and accession numbers SAMN26570417–SAMN26570460.

Contained here are files used to process the raw Illumina paired-end metagenomic and metatranscriptomic data for the above publication. These include:
* Phyloflash_workflow.sh and DNA_Phyloflash_workflow.sh - using the tool [phyloFlash](http://hrgv.github.io/phyloFlash/] to construct the SSU rRNAs (for 16S and 18S) for phylogenetic inference from the raw paired-end sequencing data.
*  DNA_phyloFlash_NTUabundance_to_phyloseq.R and TotalRNA_phyloFlash_NTUabundance_to_phyloseq.R to convert phyloFlash output data in [phyloseq](https://joey711.github.io/phyloseq/) objects in R for exoploration and visualization.
* Barplots_and_RDA_plots.R - Code for the barplots and RDA plots in the publication.

We offer these scripts in order to be as transparent as possible with our sequencing data processing. Please note that these are not intended to be run as Shell or R scripts and meant to be used as a guided walkthrough instead. Please contact me with any questions or concerns. 

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6336090.svg)](https://doi.org/10.5281/zenodo.6336090)
