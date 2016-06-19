# Methods
Genetic Map Construction
Using one unique SNP per genotype bin, we created a saturated genetic map. The genetic map was constructed using the chromosomal position of each of the SNPs as a starting point for marker ordering along the chromosomes. Each chromosome was treated as a large linkage group and each SNP was tested for linkage disequilibrium with all other SNPs with the R/QTL package (Broman et al. XXXX, Supplemental Data X) in the R statistical environment (@R-ref). The larger gaps in the map is where there is little marker information and corresponded to centromeric regions (Figure X b). These large gaps caused a small problem when ordering the markers and connecting each of the chromosomal arms in the correct order. In these ordering chromosomes X, Y, Z we used the physical position of the SNPs to connect the two arms in the correct order. The orientation of the scaffolds cannot be determined with the current methods because of we are limited by the size of the population. Therefore, scaffolds are placed between adjacent bins determined by lowest recombination probabilities.

# QTL Comparisons
To demonstrate an improvement in coverage in mapping physiological traits, we remapped two traits from @brock_floral_2010 that used the existing genetic map. As the fairest comparisons between maps, marker regression was performed using the scanone() function in RQTL with 10,000 permutations to determine the significance cutoff (Figure X).








