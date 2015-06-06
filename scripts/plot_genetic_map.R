library(qtl)
library(ggplot2)
plotMap(newmap)

plotMap
newmap
str(genetic_map)

genetic_map <- shiftmap(genetic_map)

out <- scanone(genetic_map, pheno.col = 2)
out2 <- out
head(out2)
out2$len <- 1
out2$xout   <- 1
out2$pos
str(out2)


map_plot <- ggplot(out2, aes(len, pos)) + geom_point() + facet_grid(. ~ chr)
map_plot <- map_plot + geom_errorbarh(aes(xmax = xout + 5, xmin = xout - 5, height = 0))
map_plot <- map_plot + theme(panel.background = NULL)
map_plot


map_plot <- ggplot(out2, aes(len, pos)) + facet_grid(. ~ chr)
map_plot <- map_plot + geom_segment(aes(x = len - 5, y = pos, xend = len + 5, yend = pos))
map_plot <- map_plot + theme(panel.background = NULL)
map_plot

map_plot <- ggplot(out2, aes(len, pos, ymin = 0, ymax = 40)) + facet_grid(. ~ chr)
map_plot <- map_plot + geom_segment(aes(x = len - 5, y = pos, xend = len + 5, yend = pos))
map_plot <- map_plot + scale_y_reverse()
map_plot <- map_plot + theme(panel.background = NULL, axis.text.x = element_blank(),
                             axis.line = element_blank(), panel.grid.major = element_blank(), 
                             panel.grid.minor = element_blank(), axis.ticks.x = element_blank(), 
                             axis.title.x = element_blank())
map_plot <- map_plot + ylab("Genetic Position (cM)") + ggtitle("")
map_plot


map_plot <- ggplot(out2, aes(len, pos, ymin = 0, ymax = 40)) + facet_grid(. ~ chr)
map_plot <- map_plot + geom_segment(aes(x = len - 5, y = pos, xend = len + 5, yend = pos))
map_plot <- map_plot + geom_segment(aes(x = len, y = 0, xend = len, yend = pos))
map_plot <- map_plot + scale_y_reverse()
map_plot <- map_plot + theme(panel.background = NULL, axis.text.x = element_blank(),
                             axis.line = element_blank(), panel.grid.major = element_blank(), 
                             panel.grid.minor = element_blank(), axis.ticks.x = element_blank(), 
                             axis.title.x = element_blank())
map_plot <- map_plot + ylab("Genetic Position (cM)") + ggtitle("")
map_plot

ggsave("test_genetic.pdf", map_plot)


library(qtl)
###########
# Old Map
###########
setwd("~/git.repos/brassica_meta_analysis/raw_data/")
brassica_traits <- read.cross("csvsr", genfile ="old_map_rqtl_missing_RILS_removed.csv", 
	                       phefile="Brock_2010_phenotype.csv", genotypes=c("0","2"))
head(brassica_traits)
class(brassica_traits)[1] <- "riself"
brassica_traits <- jittermap(brassica_traits)
brassica_traits
newmap <- est.map(brassica_traits,verbose=T,error.prob=.01)

pairs(jitter( as.matrix(brassica_traits$pheno) ), cex=0.6, las=1)

plot.map(brassica_traits,newmap) #some compression in this brassica_traits set
brassica_traits

brassica_traits <- replace.map(brassica_traits,newmap) #use new map
plot(brassica_traits) 
brassica_traits
plot.map(brassica_traits)

old <- scanone(brassica_traits, pheno.col = 2)
old2 <- old
head(old2)
old2$len <- 1
old2$xold   <- 1
old2$pos
str(old2)
old_map_plot <- ggplot(old2, aes(len, pos, ymin = 0, ymax = 40)) + facet_grid(. ~ chr)
old_map_plot <- old_map_plot + geom_segment(aes(x = len - 5, y = pos, xend = len + 5, yend = pos))
old_map_plot <- old_map_plot + geom_segment(aes(x = len, y = 0, xend = len, yend = pos))
old_map_plot <- old_map_plot + scale_y_reverse()
old_map_plot <- old_map_plot + theme(panel.background = NULL, axis.text.x = element_blank(),
                             axis.line = element_blank(), panel.grid.major = element_blank(), 
                             panel.grid.minor = element_blank(), axis.ticks.x = element_blank(), 
                             axis.title.x = element_blank())
old_map_plot <- old_map_plot + ylab("Genetic Position (cM)") + ggtitle("")
old_map_plot

setwd("/Users/Cody_2/git.repos/brassica_genetic_map/output")
ggsave("old_genetic_map.pdf", old_map_plot)



peak <- max(scanone_2011_flr$lod)
oldmapplotflr <- ggplot(scanone_2011_flr)
oldmapplotflr <- oldmapplotflr +  theme_bw() + geom_line(aes(x = pos, y = lod), size = 2) +
                        geom_hline(yintercept = 2.44, color = "red", size = 1) +
                        geom_segment(aes(x = pos, xend = pos), y = (peak * -0.02), yend = (peak * -0.05)) +
                        theme(text = element_text(size = 20))
                        # theme(legend.position = "none",
                        #   axis.text.x = element_text(angle = 90),
                        #   axis.line=element_line())
                        #   # panel.margin = unit(0, "cm")) +
                        # ggtitle("LOD Curves for QTLs") +
                        #xlab("Position in cM") +
                        #ylab("LOD Score") 
oldmapplotflr

###########
# New Map
###########
brassica_newmap <- read.cross("csvsr", genfile ="Brassica_F8_v1.0_gen.csv", 
	                       phefile="Brock_2010_phenotype.csv", genotypes=c("AA","BB"))
head(brassica_newmap)
class(brassica_newmap)[1] <- "riself"
brassica_newmap <- jittermap(brassica_newmap)
brassica_newmap

brassica_newmap <- est.rf(brassica_newmap)
plot.rf(brassica_newmap) 
brassica_newmap
newmap <- est.map(brassica_newmap,verbose=T,error.prob=.01)

pairs(jitter( as.matrix(brassica_newmap$pheno) ), cex=0.6, las=1)

plot.map(brassica_newmap,newmap) #some compression in this brassica_newmap set
brassica_newmap

br_phys <- read.cross("csvsr", genfile ="Brassica_F8_v2.1_gen.csv", 
                         phefile="br_blups_RQTL.csv", 
                         genotypes=c("AA","BB"), na.strings = c("-","NA"))
class(br_phys)[1] <- "riself"
br_phys
so.perm2 <- scanone(br_phys, method = "imp", n.perm = 1000) 
summary(so.perm2)
so.perm295 <- summary(so.perm2)[1] #keep 95%
# LOD thresholds (1000 permutations)
#      lod
# 5%  3.23
# 10% 2.81

?scanone
scanone_2011_STP <- scanone(br_phys, pheno.col = "X2011_STP", method = "imp", use="all.obs", chr = "A03")
plot(scanone_2011_STP, chr = "A03")

peak2 <- max(scanone_2011_STP$lod)
newplot_map <- ggplot(scanone_2011_STP)
newplot_map <- newplot_map +  theme_bw() + geom_line(aes(x = pos, y = lod), size = 2) +
                        geom_hline(yintercept = 3.08, color = "red", size = 1) +
                        geom_segment(aes(x = pos, xend = pos), y = (peak * -0.02), yend = (peak * -0.05)) +
                        theme(text = element_text(size = 20)) +
                        # theme(legend.position = "none",
                        #   axis.text.x = element_text(angle = 90),
                        #   axis.line=element_line())
                        #   # panel.margin = unit(0, "cm")) +
                        # ggtitle("LOD Curves for QTLs") +
                        xlab("Genetic Distance Along Chromosome") +
                        ylab("LOD Score") 
newplot_map


brassica_newmap <- replace.map(brassica_newmap,newmap) #use new map
plot(brassica_newmap) 
brassica_newmap
plot.map(brassica_newmap)
brassica_newmap <- shiftmap(brassica_newmap)
new <- scanone(brassica_newmap, pheno.col = 2)
new2 <- new
head(new2)
new2$len <- 1
new2$xnew   <- 1
new2$pos
str(new2)
new_map_plot <- ggplot(new2, aes(len, pos, ymin = 0, ymax = 40)) + facet_grid(. ~ chr)
new_map_plot <- new_map_plot + geom_segment(aes(x = len - 5, y = pos, xend = len + 5, yend = pos))
new_map_plot <- new_map_plot + geom_segment(aes(x = len, y = 0, xend = len, yend = pos))
new_map_plot <- new_map_plot + scale_y_reverse()
new_map_plot <- new_map_plot + theme(panel.background = NULL, axis.text.x = element_blank(),
                             axis.line = element_blank(), panel.grid.major = element_blank(), 
                             panel.grid.minor = element_blank(), axis.ticks.x = element_blank(), 
                             axis.title.x = element_blank())
new_map_plot <- new_map_plot + ylab("Genetic Position (cM)") + ggtitle("")
new_map_plot

setwd("/Users/Cody_2/git.repos/brassica_genetic_map/output")
ggsave("new_genetic_map.pdf", new_map_plot)


