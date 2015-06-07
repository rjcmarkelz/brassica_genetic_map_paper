#Brassica rapa Genetic Map v2.2

### Quickstart
Download map in RQTL format: **output/Brassica_F8_v2.2_gen.csv**

### Details
If you find a mistake or reorganize the map please issue a pull request. We are in the process of migrating all the files to open. Please be patient as some of them are tied to other publications. We are planning on making all raw files available with documentation starting the week of June 8th, 2015. 

This is a new version of the genetic map described below. This map is constructed from a total of 1527 markers from Mike Covington. These SNPs were called from RNA-seq data that was mapped to the Brassica rapa Chiffu Genome Version 1.5. Details to follow after June 8th.

###The workflow to reproduce this map:
The input file from Mike is **input/bin-genotypes_ref1.5_v0.1.1.txt** in the input directory.

I added tabs to make the file easier to parse **input/bin-genotypes_ref1.5_v0.1.1_tab.txt** also in the input directory. This is the file used as the input to the analysis pipeline I have. Run things in the following order to reproduce the map and all the figures for the genetic map portion of the paper.

#### scripts/phenotype-format.R
Format a few traits from Brock et al. (2010) for map comparison and to make RQTL objects. 

Creates: **input/brock_2010_pheno.csv**

#### scripts/snp-map-format-rqtl.R
Conversion of the genotype file into RQTL format takes **input/bin-genotypes_ref1.5_v0.1.1_tab.txt** as input. This is found in the scripts directory.

#### scripts/snp-map-construction.R
Interactive R commands to reproduce the map. Please note that many steps rely on *ripple* which has a random number component to it and may require a few calls to orderMarkers() for things to look nice. This is found in the scripts directory. I put set.seed() for where it helped to reproduce the marker order. We know the genomic location of the markers, so if chromosomes were oriented the wrong way, I reversed the marker order to match the genomic information. 

#### scripts/physical-vs-genetic.R
This compares the genetic distance to the physical distance. A few potential inversions *OR* more likely, assembly errors in the reference genome.

#### scripts/QTL-mapping-comparison.R
Genetic map plots that compare density of markers between the old map and the new map. Also, LOD score plots for two traits to compare new vs. old map QTL mapping.






















