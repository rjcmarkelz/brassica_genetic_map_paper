# Cluster bins by genotype across the RILs

```r
setwd("data")
bin.file <- "bins/bin-genotypes.scaffolds-chromosomes.2015-07-13"

bin.geno <- read.table(bin.file, header = TRUE, sep = "\t", as.is = TRUE)
bin.geno[bin.geno == "R500"]   <-  0
bin.geno[bin.geno == "IMB211"] <-  1
bin.geno[bin.geno == "HET"]    <-  NA
rownames(bin.geno) <- paste(bin.geno[, 1], bin.geno[, 2], sep = "_")

distances   <- dist(bin.geno[, 5:ncol(bin.geno)], method = "binary")
dist.m      <- as.matrix(distances)
dist.df     <- as.data.frame(dist.m)
dist.df$chr <- bin.geno$chr
dist.df$pos <- bin.geno$bin.mid

dist.df <- dist.df[!grepl('Scaffold', dist.df$chr), ]

write.table(dist.df, 'bins/bin-correlations.tsv', sep='\t', quote=F, col.names=NA)

# Can be read back in like this:
# > bin.cors <- read.table('bins/bin-correlations.tsv', sep='\t', header=T, row.names=1)
# > bin.cors$chr <- as.character(bin.cors$chr)
# > all.equal(bin.cors, dist.df)
# [1] TRUE
```
