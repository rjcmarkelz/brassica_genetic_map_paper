library(qtl)
library(ggplot2)
plotMap(newmap)

plotMap
newmap
str(genetic_map)
out <- plotMap(genetic_map)
out

out <- as.data.frame(pull.geno(genetic_map))
str(out)

out <- scanone(genetic_map, pheno.col = 2)
out2 <- head(out, 20)
out2$len <- 3

b <- ggplot(out2, aes(pos, len)) + geom_point()
b + geom_segment(aes(x = 2, y = 2, xend = 2, yend = 2))
geom_pointrange(limits)