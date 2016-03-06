#SNP Map Conversion to rqtl csvsr format
# see http://www.rqtl.org/sampledata/ for details about formats
setwd("/Users/Cody_2/git.repos/brassica_genetic_map_paper/input")
library(qtl)

# read in file and check format
snpmap <- read.delim("bin-genotypes.scaffolds-chromosomes.2015-07-13.indexed.merged",
                      header = TRUE, sep = "\t")
# snpmap2 <- read.delim("bin-genotypes_ref1.5_v0.1.1_tab.txt", header = TRUE, sep = "\t")
length(snpmap)
dim(snpmap)
# [1] 1522  132
head(snpmap)
# head(snpmap2) # double checking against other formatting
tail(snpmap)

# idx chr idx.orig chr.orig bin.mid bin.start bin.end RIL_1 RIL_103 RIL_104 RIL_113 RIL_115
# 1   1 A01        1      A01 3323536   3188282 3379775  R500    R500    R500  IMB211  IMB211
# 2   2 A01        2      A01 3444783   3423574 3469723  R500    R500    R500  IMB211  IMB211
# 3   3 A01        3      A01 3521294   3470435 3554160  R500    R500    R500  IMB211  IMB211
# 4   4 A01        4      A01 3554609   3554340 3599730  R500    R500    R500  IMB211  IMB211
# 5   5 A01        5      A01 3621755   3621755 3628191  R500    R500    R500  IMB211  IMB211
# 6   6 A01        6      A01 3679842   3642734 3707633  R500    R500    R500  IMB211  IMB211

# need to convert R500, IMB211, and HET to 2 letter genotypes as per formatting requirements
# all are imported as factors need to change to characters to replace

i <- sapply(snpmap, is.factor)
snpmap[i] <- lapply(snpmap[i], as.character)
str(snpmap)

#new data frame
snpmap_rqtl <- snpmap

snpmap_rqtl[snpmap_rqtl == "R500"] <- "AA"
snpmap_rqtl[snpmap_rqtl == "IMB211"] <- "BB"
snpmap_rqtl[snpmap_rqtl == "HET"] <- "AB"


# to add a column with the names of the markers
# make middle of bin the name of the marker
test <- snpmap_rqtl
test$id <- paste(test$chr, test$bin.mid, sep = "x")
head(test)
length(test)
colnames(test)
dim(test)
# [1] 1522  133
str(test)

# move marker name to first column and remove chr, mid, start, stop 
test_2 <- test[, c(133, 2, 5, 8:131)]

#making sure
head(test_2)
dim(test_2)
str(test_2)

#update dataframe with new rearranged columns
snpmap_rqtl <- test_2
snpmap_rqtl <- as.data.frame(snpmap_rqtl)
head(snpmap_rqtl)

#convert bp distance to Mbp
snpmap_rqtl$bin.mid <- snpmap_rqtl$bin.mid/1000000
head(snpmap_rqtl)

#rqtl csvsr format requires these two column names to be empty
names(snpmap_rqtl)
names(snpmap_rqtl)[3] <- paste("")
str(snpmap_rqtl)
snpmap_rqtl$chr <- as.character(snpmap_rqtl$chr)
names(snpmap_rqtl)[2] <- paste("")
head(snpmap_rqtl) 
# really strange----for column 3 structure(c("3.323536", "3.444783",

write.table(snpmap_rqtl, file= "snp_map_rqtl_Mbp_ref1.5.csv", row.names = FALSE,
             col.names = TRUE, sep = ",")

#end see snp-map-construction.R for next step
