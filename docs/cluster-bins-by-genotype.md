# Cluster bins by genotype across the RILs

```r
setwd("data")
bin.file <- "bins/bin-genotypes.scaffolds-chromosomes.2015-07-13"

bin.geno <- read.table(bin.file, header = TRUE, sep = "\t", as.is = TRUE)
bin.geno[bin.geno == "R500"]   <-  0
bin.geno[bin.geno == "IMB211"] <-  1
bin.geno[bin.geno == "HET"]    <-  NA

chr <- bin.geno$chr
pos <- bin.geno$bin.mid
rownames(bin.geno) <- paste(chr, pos, sep = "_")

distances   <- dist(bin.geno[, 5:ncol(bin.geno)], method = "binary")
dist.m      <- as.matrix(distances)
dist.df     <- cbind(chr, pos, as.data.frame(dist.m))
dist.df     <- dist.df[!grepl('Scaffold', dist.df$chr), ]

write.table(dist.df, 'bins/bin-correlations.tsv', sep='\t', quote=F, col.names=NA)

# Can be read back in like this:
# > bin.cors <- read.table('bins/bin-correlations.tsv', sep='\t', header=T, row.names=1)
# > bin.cors$chr <- as.character(bin.cors$chr)
# > all.equal(bin.cors, dist.df)
# [1] TRUE
```

# Plot related bins

```r
# Still have 'dist.df' defined from above and working directory set to 'data'

library(ggplot2)

out.dir <- "plots.related-bins"
dir.create(out.dir)

for (i in 3:ncol(dist.df)) {
  threshold <- quantile(dist.df[, i], 0.01)
  hits <- dist.df[, i] < threshold
  hits.dist <- dist.df[hits, i]
  hits.pos <- dist.df$pos[hits]
  hits.chr <- dist.df$chr[hits]

  marker <- names(dist.df)[i]

  if (!grepl('Scaffold', marker)) {
    chr <- sub("_\\d*", "", marker)
    chr.dist.df <- dist.df[dist.df$chr == chr, ]
    p <- ggplot(chr.dist.df) +
           geom_point(size = 1.5, aes(x = pos, y = chr.dist.df[, i],
                      color = ( chr.dist.df[, i] < threshold ))) +
           ggtitle(marker) +
           xlab("Physical Postition (Mb)") +
           ylab("Asymmetric Binary Distance") +
           guides(color = FALSE) +
           scale_colour_manual(values = c("dark grey","red")) +
           scale_x_continuous(breaks = seq(0, 40000000, 10000000),
                              labels = seq(0, 40, 10))
    ggsave(paste0(out.dir, "/", chr, ".", marker, ".png"), p, width = 11, height = 4)
  }

  if (length(unique(hits.chr)) == 1 && !grepl('Scaffold', marker)) next

  p <- ggplot(dist.df) +
         geom_point(size = 1.5, aes(x = pos, y = dist.df[, i],
                    color = ( dist.df[, i] < threshold ))) +
         facet_grid(. ~ chr, scales = "free") +
         ggtitle(marker) +
         xlab("Physical Postition (Mb)") +
         ylab("Asymmetric Binary Distance") +
         guides(color = FALSE) +
         scale_colour_manual(values = c("dark grey","red")) +
         scale_x_continuous(breaks = seq(0, 40000000, 10000000),
                            labels = seq(0, 40, 10))
  ggsave(paste0(out.dir, "/", marker, ".png"), p, width = 11, height = 4)
}
```

# Re-cluster bins after they have been rearranged to assign scaffolds and fix misassemblies

```r
setwd("data")
bin.file <- "bins/bin-genotypes.scaffolds-chromosomes.2015-07-13.indexed"

bin.geno <- read.table(bin.file, header = TRUE, sep = "\t", as.is = TRUE)
bin.geno[bin.geno == "R500"]   <-  0
bin.geno[bin.geno == "IMB211"] <-  1
bin.geno[bin.geno == "HET"]    <-  NA

chr <- bin.geno$chr
idx <- bin.geno$idx
rownames(bin.geno) <- paste(chr, idx, sep = "_")

distances   <- dist(bin.geno[, 8:(ncol(bin.geno) - 1)], method = "binary")
dist.m      <- as.matrix(distances)
dist.df     <- cbind(chr, idx, as.data.frame(dist.m))

write.table(dist.df, 'bins/bin-correlations.rearranged.tsv', sep='\t', quote=F, col.names=NA)
```

# Plot pairwise distances for bins after they have been rearranged to assign scaffolds and fix misassemblies

```r
library(ggplot2)
library(grid)
library(reshape2)

dist.melt <- melt(dist.df, id.vars = c('chr', 'idx'))
dist.melt <- cbind(dist.melt, colsplit(dist.melt$variable, '_', c('chr2', 'idx2')))
dist.melt$chr2 <- factor(dist.melt$chr2, levels = rev(levels(dist.melt$chr)))

p <- ggplot(dist.melt, aes(x = idx, y = idx2)) +
       geom_tile(aes(fill = value)) +
       scale_fill_gradientn(limits = c(0, 1), colours = c('dodgerblue4', 'dodgerblue', 'aliceblue', 'white')) +
       facet_grid(chr2 ~ chr, scales = 'free', space = 'free') +
       theme_bw() +
       theme(
         axis.text = element_blank(),
         axis.ticks = element_blank(),
         axis.title = element_blank(),
         panel.border = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         panel.margin = unit(0.1, "lines"),
         legend.background = element_rect(colour = "darkgray", fill='white'),
         legend.box='horizontal',
         legend.direction='vertical',
         legend.justification = c(0, 1),
         legend.position = c(0.025,0.95)
       ) +
       guides(
         fill = guide_colorbar(barwidth = 2,
                               reverse = TRUE,
                               title = 'Asymmetric\nBinary\nDistance\n',
                               title.hjust = 0.5,
                               title.position = 'top')
       ) +
       scale_x_continuous(expand = c(0, 0)) +
       scale_y_continuous(expand = c(0, 0))

ggsave('bins/bin-distances.rearranged.png', p, width = 8, height = 8)
```
