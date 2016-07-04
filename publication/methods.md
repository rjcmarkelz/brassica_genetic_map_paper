# Methods
## Genotyping 
Individuals in the population were genotyped first using the reference set of sSNPs that were called as part of Devissetty et al. 2014 pipeline. However, once merging all the outputs for the RILs there were inconsistencies suggesting that the parents of the population had more then one parent. The R500 SNPs were consistent across all the individuals whereas the IMB211 SNPs were not suggesting that the parental seed stock has cross pollination contamination. There was not any evidence that the R500 seed stock was contaminated. This led us to approach this using an alternative genotyping strategy by combining all the individual replicates per RIL into a single file to call SNPs compared to the reference genome. Each SNP was then filtered by comparing it to the genotype of the R500 parent. At each SNP, each RIL was genotyped as R500 or the alternate allele. All these SNPs needed to match up with the one another. At each step there were quality scores and filters that were applied using custom Perl scripts (found here). Only SNPs that met our quality scores were considered. 

## Genotypic Bin Creation 
All of the SNPs were then assembled along the chromosomes according to their genomic location. Aligning all the individuals in the population allowed us to calculate unique recombinations that are present in the population. These unique bins were determined by custom Perl scripts (found here). SNPs located in the middle of the genotypic bin were selected to be used for placement of scaffolds and the creation of the genetic map. 

## Genetic Map Construction
Using one unique SNP per genotype bin, we created a saturated genetic map. The genetic map was constructed using the chromosomal position of each of the SNPs as a starting point for marker ordering along the chromosomes. Each chromosome was treated as a large linkage group and each SNP was tested for linkage disequilibrium with all other SNPs with the R/QTL package (Broman et al. XXXX, Supplemental Data X) in the R statistical environment (@R-ref). The larger gaps in the map is where there is little marker information and corresponded to centromeric regions (Figure X b). These large gaps caused a small problem when ordering the markers and connecting each of the chromosomal arms in the correct order. In these ordering chromosomes X, Y, Z we used the physical position of the SNPs to connect the two arms in the correct order. The orientation of the scaffolds cannot be determined with the current methods because of we are limited by the size of the population. Therefore, scaffolds are placed between adjacent bins determined by lowest recombination probabilities.

## QTL Comparisons
To demonstrate an improvement in coverage in mapping physiological traits, we remapped two traits from @brock_floral_2010 that used the existing genetic map. As the fairest comparisons between maps, marker regression was performed using the scanone() function in RQTL with 10,000 permutations to determine the significance cutoff (Figure X).







