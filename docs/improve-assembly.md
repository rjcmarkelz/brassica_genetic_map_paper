
Make a version of the combined scaffolds and chromosomes bin genotype file with columns for original chromosome and order on that chromosome for each bin prior to making any rearrangements to improve the assembly.

```r
setwd('data/bin/')
bin.genotypes <- read.table('bin-genotypes.scaffolds-chromosomes.2015-07-13', header=T)

chr.orig <- NaN
idx.orig <- NaN
bin.genotypes <- cbind(chr.orig, idx.orig, bin.genotypes)

for (chr in unique(bin.genotypes$chr)) {
  if (grepl('Scaffold', chr)) next
  chr.current <- bin.genotypes$chr == chr
  bin.genotypes$chr.orig[chr.current] <- as.character(bin.genotypes$chr[chr.current])
  bin.genotypes$idx.orig[chr.current] <- seq(1, length(bin.genotypes$idx.orig[chr.current]))
}
write.table(bin.genotypes, 'bin-genotypes.scaffolds-chromosomes.2015-07-13.indexed', quote=F, sep='\t', row.names=F)
```
