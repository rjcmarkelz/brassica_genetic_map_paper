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
out2 <- out
head(out2)
out2$len <- 1
out2$xout   <- 1
out2$pos
str(out2)


b <- ggplot(out2, aes(len, pos)) + geom_point() + facet_grid(. ~ chr)
b + geom_errorbarh(aes(xmax = xout + 5, xmin = xout - 5, height = 0))







