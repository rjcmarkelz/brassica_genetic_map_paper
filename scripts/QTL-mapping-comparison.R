
###########
# Old Map
###########
library(qtl)
library(ggplot2)
setwd("/Users/Cody_2/git.repos/brassica_genetic_map_paper/input")
brassica_traits <- read.cross("csvsr", genfile ="old_map_rqtl.csv", 
	                       phefile="Brock_2010_pheno.csv", genotypes=c("0","2"))
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

#change the names of the chromosomes
old_map <- pull.map(brassica_traits)
str(brassica_traits)
names(brassica_traits$geno) <- paste(c("A01", "A02", "A03", "A04", "A05", "A06", "A07", "A08", "A09", "A10"))
names(brassica_traits$geno)

old <- scanone(brassica_traits, pheno.col = 1, method = "imp")
old
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
old_map_plot <- old_map_plot + theme(panel.background = element_blank(), axis.text.x = element_blank(),
                             axis.line = element_blank(), panel.grid.major = element_blank(), 
                             panel.grid.minor = element_blank(), axis.ticks.x = element_blank(), 
                             axis.title.x = element_blank())
old_map_plot <- old_map_plot + ylab("Genetic Position (cM)") + ggtitle("")
old_map_plot

setwd("/Users/Cody_2/git.repos/brassica_genetic_map/output")
ggsave("old_genetic_map.pdf", old_map_plot)


oldleaf <- scanone(brassica_traits, pheno.col = 1, method = "imp", chr = "A06")
plot(oldleaf)

peak <- 7
oldmapplotleaf <- ggplot(oldleaf)
oldmapplotleaf <- oldmapplotleaf +  theme_bw() + scale_y_continuous(limits=c(0, 7.5)) + 
                        geom_line(aes(x = pos, y = lod), size = 2) +
                        geom_hline(yintercept = 2.44, color = "red", size = 1) +
                        geom_segment(aes(x = pos, xend = pos), y = (peak * -0.02), yend = (peak * -0.05)) +
                        theme(text = element_text(size = 20)) +
                        # theme(legend.position = "none",
                        #   axis.text.x = element_text(angle = 90),
                        #   axis.line=element_line())
                        #   # panel.margin = unit(0, "cm")) +
                        # ggtitle("LOD Curves for QTLs") +
                        xlab("Genetic Position (cM)") +
                        ylab("LOD Score") 
oldmapplotleaf

oldmapplotflr <- scanone(brassica_traits, pheno.col = 2, method = "imp", chr = "A10")
plot(oldmapplotflr)
peak2 <- 9
oldplot_flr <- ggplot(oldmapplotflr)
oldplot_flr <- oldplot_flr +  theme_bw() + scale_y_continuous(limits=c(0, 13)) + 
                        geom_line(aes(x = pos, y = lod), size = 2) +
                        geom_hline(yintercept = 3.08, color = "red", size = 1) +
                        geom_segment(aes(x = pos, xend = pos), y = (peak * -0.02), yend = (peak * -0.05)) +
                        theme(text = element_text(size = 20)) +
                        # theme(legend.position = "none",
                        #   axis.text.x = element_text(angle = 90),
                        #   axis.line=element_line())
                        #   # panel.margin = unit(0, "cm")) +
                        # ggtitle("LOD Curves for QTLs") +
                        xlab("Genetic Distance (cM)") +
                        ylab("LOD Score") 
oldplot_flr


###########
# New Map
###########
setwd("/Users/Cody_2/git.repos/brassica_genetic_map/output")
brassica_newmap <- read.cross("csvsr", genfile ="snp_map_rqtl_Mbp_ref1.5_cross_output_gen.csv", 
                           phefile="snp_map_rqtl_Mbp_ref1.5_cross_output_phe.csv", 
                           genotypes=c("AA","BB"), 
                           na.strings=c("NA","-"))
head(brassica_newmap)
class(brassica_newmap)[1] <- "riself"
brassica_newmap <- jittermap(brassica_newmap)
brassica_newmap

newmapplotleaf_perm <- scanone(brassica_newmap, pheno.col = 1, method = "imp") 
summary(so.perm2)
# so.perm295 <- summary(so.perm2)[1] #keep 95%
# # LOD thresholds (1000 permutations)
# #      lod
# # 5%  3.23
# # 10% 2.81

newmapplotleaf <- scanone(brassica_newmap, pheno.col = 1, method = "imp", chr = "A06")
plot(newmapplotleaf)

peak2 <- 7
newplot_map <- ggplot(newmapplotleaf)
newplot_map <- newplot_map +  theme_bw() + scale_y_continuous(limits=c(0, 7.5)) + 
                        geom_line(aes(x = pos, y = lod), size = 2) +
                        geom_hline(yintercept = 3.08, color = "red", size = 1) +
                        geom_segment(aes(x = pos, xend = pos), y = (peak * -0.02), yend = (peak * -0.05)) +
                        theme(text = element_text(size = 20)) +
                        # theme(legend.position = "none",
                        #   axis.text.x = element_text(angle = 90),
                        #   axis.line=element_line())
                        #   # panel.margin = unit(0, "cm")) +
                        # ggtitle("LOD Curves for QTLs") +
                        xlab("Genetic Distance (cM)") +
                        ylab("LOD Score") 
newplot_map


newmapplotflr <- scanone(brassica_newmap, pheno.col = 2, method = "imp", chr = "A10")
plot(newmapplotflr)
peak2 <- 9
newplot_flr <- ggplot(newmapplotflr)
newplot_flr <- newplot_flr +  theme_bw() + scale_y_continuous(limits=c(0, 13)) + 
                        geom_line(aes(x = pos, y = lod), size = 2) +
                        geom_hline(yintercept = 3.08, color = "red", size = 1) +
                        geom_segment(aes(x = pos, xend = pos), y = (peak * -0.02), yend = (peak * -0.05)) +
                        theme(text = element_text(size = 20)) +
                        # theme(legend.position = "none",
                        #   axis.text.x = element_text(angle = 90),
                        #   axis.line=element_line())
                        #   # panel.margin = unit(0, "cm")) +
                        # ggtitle("LOD Curves for QTLs") +
                        xlab("Genetic Distance (cM)") +
                        ylab("LOD Score") 
newplot_flr


new <- scanone(brassica_newmap, pheno.col = 2)
new2 <- new
head(new2)
new2$len <- 1
new2$xnew <- 1
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


