
Make a version of the combined scaffolds and chromosomes bin genotype file with columns for original chromosome and order on that chromosome for each bin prior to making any rearrangements to improve the assembly.

```r
setwd('data/bins/')
bin.genotypes <- read.table('bin-genotypes.scaffolds-chromosomes.2015-07-13', header=T)

chr.new <- 'NA'
idx.orig <- NaN
comments <- ''
bin.genotypes <- cbind(chr.new, idx.orig, bin.genotypes, comments)
bin.genotypes$chr.new <- as.character(bin.genotypes$chr.new)
bin.genotypes$comments <- as.character(bin.genotypes$comments)

for (chr in unique(bin.genotypes$chr)) {
  if (grepl('Scaffold', chr)) next
  chr.current <- bin.genotypes$chr == chr
  bin.genotypes$chr.new[chr.current] <- as.character(bin.genotypes$chr[chr.current])
  bin.genotypes$idx.orig[chr.current] <- seq(1, length(bin.genotypes$idx.orig[chr.current]))
}
write.table(bin.genotypes, 'bin-genotypes.scaffolds-chromosomes.2015-07-13.indexed', quote=F, sep='\t', row.names=F)
```

Function to return the top 10 similar bins:

```r
similar.bins <- function(bin, df=dist.df, top.n=10) {
  bin.df <- df[, bin, drop=F]
  bin.df <- bin.df[order(bin.df[, 1]), , drop=F]
  return(bin.df[1:top.n, , drop=F])
}
```

Usage example:

```r
similar.bins('Scaffold000100_416921')
```

Output example:

                 Scaffold000100_416921
    A02_21152172            0.00000000
    A06_12814512            0.00000000
    A06_11042216            0.01333333
    A06_14144472            0.01351351
    A06_10175707            0.02666667
    A06_14379326            0.03947368
    A06_9978506             0.04000000
    A06_14927290            0.05263158
    A06_15466097            0.06578947
    A06_9706189             0.06666667

After rearranging bins to place Scaffolds and fix genome misassemblies, add a column with the new chromosome index (`idx.new`):

```r
setwd('data/bins/')
bin.genotypes <- read.table('bin-genotypes.scaffolds-chromosomes.2015-07-13.indexed', header=T, fill=T)

idx.new <- NaN
bin.genotypes <- cbind(idx.new, bin.genotypes)

for (chr in unique(bin.genotypes$chr.new)) {
  chr.current <- which(bin.genotypes$chr.new == chr)
  bin.genotypes$idx.new[chr.current] <- seq(1, length(chr.current))
}
write.table(bin.genotypes, 'bin-genotypes.scaffolds-chromosomes.2015-07-13.indexed', quote=F, sep='\t', row.names=F)
```
