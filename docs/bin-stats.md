


```r
setwd("data")
bin.file <- "bins/bin-genotypes.scaffolds-chromosomes.2015-07-13.indexed"

bin.geno <- read.table(bin.file, header = TRUE, sep = "\t", as.is = TRUE)

flagged <- rep(FALSE, nrow(bin.geno))
flagged[bin.geno$chr != bin.geno$chr.orig] <- TRUE

for (i in 1:(nrow(bin.geno) - 1)) {
  if (!is.na(bin.geno$idx.orig[i]) && !is.na(bin.geno$idx.orig[i + 1]) &&
      bin.geno$idx.orig[i] + 1 != bin.geno$idx.orig[i + 1] &&
      bin.geno$chr[i] == bin.geno$chr[i + 1]) {
    flagged[max(1, (i - 4)):min((i + 5), nrow(bin.geno))] <- TRUE
  }
}

write.table(bin.geno[flagged, 1:7], 'bins/flagged.genotypes.tsv', sep='\t', quote=F, col.names=NA)
```











