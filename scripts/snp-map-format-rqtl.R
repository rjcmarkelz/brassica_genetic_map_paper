#SNP Map Conversion to rqtl csvsr format
# see http://www.rqtl.org/sampledata/ for details about formats

setwd("/Users/Cody_2/git.repos/brassica_genetic_map_paper/input")

snpmap <- read.delim("bin-genotypes_ref1.5_v0.1.1_tab.txt", header = TRUE, sep = "\t")
length(snpmap)
dim(snpmap)
# [1] 1527  128
head(snpmap)
tail(snpmap)

# example
#     chr  bin.mid bin.start  bin.end  RIL_1 RIL_103 RIL_104 RIL_113 RIL_115 RIL_12 RIL_123
# 1522 A10 16080749  16080747 16082372 IMB211  IMB211  IMB211  IMB211  IMB211   R500  IMB211
# 1523 A10 16115802  16084408 16156351 IMB211  IMB211  IMB211  IMB211  IMB211   R500  IMB211
# 1524 A10 16192145  16167764 16208022 IMB211  IMB211  IMB211  IMB211  IMB211   R500  IMB211


#need to convert R500, IMB211, and HET to 2 letter genotypes as per formatting requirements
#all are imported as factors need to change to characters to replace
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
# [1] 1527  129
str(test)

# move marker name to first column and remove chr, mid, start, stop 
test_2 <- test[,c(129, 1, 2, 5:128)]
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
names(snpmap_rqtl)[3] <- paste("")
names(snpmap_rqtl)[2] <- paste("")
head(snpmap_rqtl)

write.table(snpmap_rqtl, file= "snp_map_rqtl_Mbp_ref1.5.csv", row.names = FALSE, col.names = TRUE, sep = ",")
#end  
