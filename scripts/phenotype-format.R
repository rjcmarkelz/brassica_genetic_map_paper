# Phenotype Data Conversion to rqtl csvsr format
# see http://www.rqtl.org/sampledata/ for details about formats

library(qtl)
library(readxl)

# FileS1 From here:
# http://www.genetics.org/content/186/4/1451/rel-suppl/698981c0c2877a54/suppl/DC1
setwd("/Users/Cody_2/git.repos/brassica_genetic_map_paper/input")

#read in blup data
brock_blups <- read_excel("FileS1.xlsx", sheet = "Greenhouse07", na = "NA")
head(brock_blups)
str(brock_blups)

# name things the same as the genetic map, replace Line column with RIL
brock_blups$RILs <- sub("(Ind-)(\\d+)", "RIL_\\2", brock_blups$Line)
dim(brock_blups)

# keep two traits for remapping comparison
brock_blups <- brock_blups[,c(3,6,14)]
head(brock_blups)

# reformat for RQTL, csvsr format to pull in same phenotype data 
# with different genetic maps for comparison
brock_t <- as.data.frame(t(brock_blups))
row.names(brock_t)[3] <- paste("id")

write.table(brock_t, "brock_2010_pheno.csv", row.names = T, col.names = F, sep = ",")

# note there will be some genotypes that are not represented with the new map that will be 
# removed in subsequent analyses with the dense genetic map