#Materials and Methods

###Genetic Map Construction

Using one unique SNP per genotype bin, we created a saturated genetic map.

The genetic map was constructed using the chromosomal position of each of the SNPs as a starting point for marker ordering along the chromosomes.

The genetic map construction was done using the R statistical environment with the R/QTL package (Broman et al. XXXX, Broman genetic map paper, Supplemental Data X).

In brief, each chromosome was treated as a large linkage group and each SNP was tested for linkage disequilibrium with all other SNPs to find and move misplaced markers.

The larger gaps in the map are where we had little marker information and corresponded to centromeric regions.

These large gaps caused a small problem when ordering the markers and connecting each of the chromosomal arms in the correct order.

On chromosomes X, Y, Z we used the physical position of the SNPs to connect the two arms in the correct order. 

###QTL Mapping

To demonstrate an improvement in coverage in mapping physiological traits, we remapped three traits from Brock et al. 2010 that used the old genetic map without physical positions of the markers.

This was performed using the scanone() function in RQTL with 10,000 permutations to determine the significance cutoff (see Supplemental Data X.)

Comparisons of the old vs. new map were visualized using ggplot2 (Wickham XXXX).






