library(qtl)
setwd("/Users/Cody_2/git.repos/brassica_genetic_map/input")

brassica_genes <- read.cross("csvsr", genfile ="snp_map_rqtl_Mbp_ref1.5.csv", 
	                       phefile="phenotype.csv", 
	                       genotypes=c("AA","BB"), 
	                       na.strings=c("NA","AB"))

head(brassica_genes)
summary(brassica_genes)

# change to RIL population
class(brassica_genes)[1] <- "riself"
brassica_genes <- jittermap(brassica_genes)
brassica_genes

# takes about a minute to estimate the map
newmap <- est.map(brassica_genes, error.prob= 0.0005)
summaryMap(newmap)
#         n.mar  length ave.spacing max.spacing
# A01       175  7171.3        41.2       966.8
# A02       125  3169.3        25.6       966.8
# A03       216   135.6         0.6         4.6
# A04       106    63.2         0.6         2.1
# A05        91  2384.6        26.5       966.8
# A06       168    93.9         0.6         2.9
# A07       166   489.5         3.0       139.9
# A08        84  2127.3        25.6       966.8
# A09       224   175.5         0.8        14.4
# A10       172   121.4         0.7        21.2
# overall  1527 15931.6        10.5       966.8

plotMap(newmap) # saved output _1
plot.rf(brassica_genes) #saved output _1

#find duplicates and drop them
dup <- findDupMarkers(brassica_genes, exact.only=FALSE)
dup
head(dup)
length(dup)
str(dup)
# 57

totmar(brassica_genes)
# 1527
brassica_genes <- drop.markers(brassica_genes, unlist(dup))
totmar(brassica_genes)
# 1465

#identify markers with possible segregation distortion
gt <- geno.table(brassica_genes)
seg_dist <- gt[gt$P.value < 0.05/totmar(brassica_genes),]
seg_dist
dim(seg_dist)
# 33 all on chromosome 3, will not drop these, but make note
#              chr missing AA  BB      P.value
# A03x19849963 A03       0 38  86 1.628658e-05
# A03x19928814 A03       0 36  88 3.015843e-06
# A03x20027810 A03       0 37  87 7.117887e-06
# A03x20139956 A03       0 36  88 3.015843e-06 ...... truncated output

# reestimate map with dropped duplicate markers
newmap <- est.map(brassica_genes, error.prob= 0.0005)
summaryMap(newmap)
plotMap(newmap) # saved _2
brassica_genes <- replace.map(brassica_genes,newmap)
plot.rf(brassica_genes) # saved _2



# now to get into dirty marker dropping mode
# make a copy so I can go back and not have to reestimate the map
brassica_genes_reduced_2 <- brassica_genes
rffull <- pull.rf(brassica_genes_reduced_2)

#############
# warm up with A08
#############
plot.rf(brassica_genes_reduced_2, chr='A08')
rf8 <- pull.rf(brassica_genes_reduced_2, chr = 'A08')
chr8 <- markernames(brassica_genes_reduced_2, chr='A08')
chr8
#plotting is the easiest

plot(rf8, chr8[7], bandcol="gray70", ylim=c(0,1), alternate.chrid=TRUE)
plot(rf8, chr8[8], bandcol="gray70", ylim=c(0,1), alternate.chrid=TRUE)
#problem markers are 7 and 8
plot(rf8, chr8[9], bandcol="gray70", ylim=c(0,1), alternate.chrid=TRUE)
plot(rf8, chr8[30], bandcol="gray70", ylim=c(0,1), alternate.chrid=TRUE)

#remove 7 and 8, could also belong on chromosome 9 see plots below
plot(rffull, chr8[7], bandcol="gray70", ylim=c(0,1), alternate.chrid=TRUE)
plot(rffull, chr8[8], bandcol="gray70", ylim=c(0,1), alternate.chrid=TRUE)
brassica_genes_reduced_2 <- drop.markers(brassica_genes_reduced_2,
	c("A08x4875175",  "A08x4952134"))


brassica_genes_reduced_2 <- orderMarkers(brassica_genes_reduced_2, chr=c('A08'), 
	                        window=5, use.ripple=TRUE, maxit=4000, 
	                        error.prob=0.0001)
plotMap(brassica_genes_reduced_2, chr='A08') 
plot.rf(brassica_genes_reduced_2, chr='A08')
rf8 <- pull.rf(brassica_genes_reduced_2, chr = 'A08')
chr8 <- markernames(brassica_genes_reduced_2, chr='A08')
chr8
#still something funny
plot(rffull, chr8[1], bandcol="gray70", ylim=c(0,1), alternate.chrid=TRUE)
plot(rffull, chr8[2], bandcol="gray70", ylim=c(0,1), alternate.chrid=TRUE)
plot(rffull, chr8[3], bandcol="gray70", ylim=c(0,1), alternate.chrid=TRUE)
# markers "A08x1515279"  "A08x1485026" belong on chromosome 2
# remove for now
brassica_genes_reduced_2 <- drop.markers(brassica_genes_reduced_2,
	c("A08x1515279",  "A08x1485026"))

brassica_genes_reduced_2 <- orderMarkers(brassica_genes_reduced_2, chr=c('A08'), 
	                        window=5, use.ripple=TRUE, maxit=4000, 
	                        error.prob=0.0001)
plotMap(brassica_genes_reduced_2, chr='A08') 
plot.rf(brassica_genes_reduced_2, chr='A08')


plotMap(brassica_genes_reduced_2)

#############
# Onto A05
#############
plot.rf(brassica_genes_reduced_2, chr='A05')
rf5 <- pull.rf(brassica_genes_reduced_2, chr = 'A05')
chr5 <- markernames(brassica_genes_reduced_2, chr='A05')
chr5

plot(rf5, chr5[22], bandcol="gray70", ylim=c(0,1), alternate.chrid=TRUE)
plot(rf5, chr5[23], bandcol="gray70", ylim=c(0,1), alternate.chrid=TRUE)
plot(rf5, chr5[24], bandcol="gray70", ylim=c(0,1), alternate.chrid=TRUE)
plot(rf5, chr5[25], bandcol="gray70", ylim=c(0,1), alternate.chrid=TRUE)
plot(rf5, chr5[26], bandcol="gray70", ylim=c(0,1), alternate.chrid=TRUE)
plot(rf5, chr5[27], bandcol="gray70", ylim=c(0,1), alternate.chrid=TRUE)
plot(rf5, chr5[28], bandcol="gray70", ylim=c(0,1), alternate.chrid=TRUE)

plot(rffull, chr5[22], bandcol="gray70", ylim=c(0,1), alternate.chrid=TRUE)
plot(rffull, chr5[23], bandcol="gray70", ylim=c(0,1), alternate.chrid=TRUE)

plot(rffull, chr5[24], bandcol="gray70", ylim=c(0,1), alternate.chrid=TRUE)
plot(rffull, chr5[25], bandcol="gray70", ylim=c(0,1), alternate.chrid=TRUE)
plot(rffull, chr5[26], bandcol="gray70", ylim=c(0,1), alternate.chrid=TRUE)
plot(rffull, chr5[27], bandcol="gray70", ylim=c(0,1), alternate.chrid=TRUE)

# this batch of bad markers belongs on chr A08
# "A05x8149908" "A05x8245839"  "A05x8277002"  "A05x8327143"
brassica_genes_reduced_2 <- drop.markers(brassica_genes_reduced_2,
	c("A05x8149908", "A05x8245839",  "A05x8277002",  "A05x8327143"))

brassica_genes_reduced_2 <- orderMarkers(brassica_genes_reduced_2, chr=c('A05'), 
	                        window=5, use.ripple=TRUE, maxit=4000, 
	                        error.prob=0.0001)
plotMap(brassica_genes_reduced_2, chr='A05') 
plot.rf(brassica_genes_reduced_2, chr='A05')


plotMap(brassica_genes_reduced_2)

#############
# A02
#############
plot.rf(brassica_genes_reduced_2)
plot.rf(brassica_genes_reduced_2, chr='A02')
rf2 <- pull.rf(brassica_genes_reduced_2, chr = 'A02')
chr2 <- markernames(brassica_genes_reduced_2, chr='A02')
chr2
# marker "A02x21152172" is strange and belongs on A06
plot(rf2, chr2[76], bandcol="gray70", ylim=c(0,1), alternate.chrid=TRUE)
plot(rf2, chr2[77], bandcol="gray70", ylim=c(0,1), alternate.chrid=TRUE)
plot(rf2, chr2[78], bandcol="gray70", ylim=c(0,1), alternate.chrid=TRUE)
plot(rf2, chr2[79], bandcol="gray70", ylim=c(0,1), alternate.chrid=TRUE)

plot(rffull, chr2[78], bandcol="gray70", ylim=c(0,1), alternate.chrid=TRUE)

# marker "A02x23373611" is strange and belongs on chromosome A10
plot(rf2, chr2[91], bandcol="gray70", ylim=c(0,1), alternate.chrid=TRUE)
plot(rf2, chr2[92], bandcol="gray70", ylim=c(0,1), alternate.chrid=TRUE)
plot(rf2, chr2[93], bandcol="gray70", ylim=c(0,1), alternate.chrid=TRUE)

plot(rffull, chr2[93], bandcol="gray70", ylim=c(0,1), alternate.chrid=TRUE)

brassica_genes_reduced_2 <- drop.markers(brassica_genes_reduced_2,
	c("A02x21152172", "A02x23373611"))

# may need to come back to 2, it is still really large ~200 cM
brassica_genes_reduced_2 <- orderMarkers(brassica_genes_reduced_2, chr=c('A02'), 
	                        window=5, use.ripple=TRUE, maxit=4000, 
	                        error.prob=0.0001)
plotMap(brassica_genes_reduced_2, chr='A02') 
plot.rf(brassica_genes_reduced_2, chr='A02')

plotMap(brassica_genes_reduced_2)

#############
# A01
#############
plot.rf(brassica_genes_reduced_2, chr='A01')
rf1 <- pull.rf(brassica_genes_reduced_2, chr = 'A01')
chr1 <- markernames(brassica_genes_reduced_2, chr='A01')
chr1

# marker "A01x8900522" belongs on A09
plot(rf1, chr1[58], bandcol="gray70", ylim=c(0,1), alternate.chrid=TRUE)
plot(rf1, chr1[59], bandcol="gray70", ylim=c(0,1), alternate.chrid=TRUE)
plot(rffull, chr1[58], bandcol="gray70", ylim=c(0,1), alternate.chrid=TRUE)

# marker "A01x10301489" belongs on A07
plot(rf1, chr1[79], bandcol="gray70", ylim=c(0,1), alternate.chrid=TRUE)
plot(rffull, chr1[79], bandcol="gray70", ylim=c(0,1), alternate.chrid=TRUE)

# marker "A01x11450244" belongs on A04
plot(rf1, chr1[84], bandcol="gray70", ylim=c(0,1), alternate.chrid=TRUE)
plot(rffull, chr1[84], bandcol="gray70", ylim=c(0,1), alternate.chrid=TRUE)

# marker "A01x16479326" belongs on A02
plot(rf1, chr1[92], bandcol="gray70", ylim=c(0,1), alternate.chrid=TRUE)
plot(rffull, chr1[92], bandcol="gray70", ylim=c(0,1), alternate.chrid=TRUE)

brassica_genes_reduced_2 <- drop.markers(brassica_genes_reduced_2,
	c("A01x8900522", "A01x10301489", "A01x11450244", "A01x16479326"))

#change window size to 2 and it falls into place nicely
brassica_genes_reduced_2 <- orderMarkers(brassica_genes_reduced_2, chr=c('A01'), 
	                        window=2, use.ripple=TRUE, maxit=4000, 
	                        error.prob=0.0001)
plotMap(brassica_genes_reduced_2, chr='A01') 
plot.rf(brassica_genes_reduced_2, chr='A01')



# save output progress
setwd("/Users/Cody_2/git.repos/brassica_genetic_map/Output")
write.cross(brassica_genes_reduced_2, format= "csvsr", 
	         filestem="snp_map_rqtl_Mbp_ref1.5_cross_output")

plotMap(brassica_genes_reduced_2) 
plot.rf(brassica_genes_reduced_2)


#############
# A09 ripple marker order
#############
plotMap(brassica_genes_reduced_2, chr='A09') 
plot.rf(brassica_genes_reduced_2, chr = 'A09')

?orderMarkers
brassica_genes_reduced_2 <- orderMarkers(brassica_genes_reduced_2, chr = c('A09'), 
	                        window = 2, use.ripple = TRUE, maxit = 4000, 
	                        error.prob = 0.0001)

plotMap(brassica_genes_reduced_2, chr='A09') 
plot.rf(brassica_genes_reduced_2, chr = 'A09')


#############
# A10 ripple marker order
#############
plotMap(brassica_genes_reduced_2, chr='A10') 
plot.rf(brassica_genes_reduced_2, chr = 'A10')

?orderMarkers
brassica_genes_reduced_2 <- orderMarkers(brassica_genes_reduced_2, chr = c('A10'), 
	                        window = 2, use.ripple = TRUE, maxit = 4000, 
	                        error.prob = 0.0001)

plotMap(brassica_genes_reduced_2, chr='A10') 
plot.rf(brassica_genes_reduced_2, chr = 'A10')

# save output progress
setwd("/Users/Cody_2/git.repos/brassica_genetic_map/Output")
write.cross(brassica_genes_reduced_2, format= "csvsr", 
	         filestem="snp_map_rqtl_Mbp_ref1.5_cross_output")

##########################
##########################
# Fine tuning
##########################
##########################
summaryMap(brassica_genes_reduced_2)
#         n.mar length ave.spacing max.spacing
# A01       163   94.0         0.6         3.1
# A02       119  121.3         1.0        36.6
# A03       207  135.6         0.7         4.6
# A04       100   63.3         0.6         2.1
# A05        84   91.1         1.1        26.3
# A06       158   93.9         0.6         2.9
# A07       160  106.4         0.7         3.1
# A08        77   80.7         1.1        29.1
# A09       216  154.7         0.7        14.5
# A10       167  104.1         0.6         2.5
# overall  1451 1045.0         0.7        36.6

#############
# Back to A02
#############
plotMap(brassica_genes_reduced_2)
plot.rf(brassica_genes_reduced_2, chr='A02')

rffull <- pull.rf(brassica_genes_reduced_2)
rf2 <- pull.rf(brassica_genes_reduced_2, chr = 'A02')
chr2 <- markernames(brassica_genes_reduced_2, chr='A02')
chr2

plot(rf2, chr2[20], bandcol="gray70", ylim=c(0,1), alternate.chrid=TRUE)
plot(rffull, chr2[70], bandcol="gray70", ylim=c(0,1), alternate.chrid=TRUE)

# just needed to reorder the markers on the top of A02
brassica_genes_reduced_2 <- orderMarkers(brassica_genes_reduced_2, chr=c('A02'), 
	                        window=2,  maxit=4000, 
	                        error.prob=0.0001)
plotMap(brassica_genes_reduced_2, chr='A02') 
plot.rf(brassica_genes_reduced_2, chr='A02')


#############
# Back to A08
#############
plot.rf(brassica_genes_reduced_2, chr='A08')
rf8 <- pull.rf(brassica_genes_reduced_2, chr = 'A08')
chr8 <- markernames(brassica_genes_reduced_2, chr='A08')
chr8
# need to flip order of markers
brassica_genes_reduced_2 <- orderMarkers(brassica_genes_reduced_2, chr=c('A08'), 
	                        window=2, use.ripple = TRUE, maxit=4000, 
	                        error.prob=0.0001)
plotMap(brassica_genes_reduced_2, chr='A08') 
plot.rf(brassica_genes_reduced_2, chr='A08')


#############
# Back to A05
#############
rf5 <- pull.rf(brassica_genes_reduced_2, chr = 'A05')
chr5 <- markernames(brassica_genes_reduced_2, chr='A05')
chr5
plot.rf(brassica_genes_reduced_2, chr='A05')
plotMap(brassica_genes_reduced_2, chr='A05') 
brassica_genes_reduced_2 <- orderMarkers(brassica_genes_reduced_2, chr=c('A05'), 
	                        window=2, use.ripple = TRUE, maxit=4000, 
	                        error.prob=0.0001)
plotMap(brassica_genes_reduced_2, chr='A05') 
plot.rf(brassica_genes_reduced_2, chr='A05')

#cannot make any improvements on A05

# final plots
plotMap(brassica_genes_reduced_2) 
plot.rf(brassica_genes_reduced_2)

# do the marker orders make sense?
chr1 <- markernames(brassica_genes_reduced_2, chr='A01')
chr1
chr2 <- markernames(brassica_genes_reduced_2, chr='A02')
chr2
chr3 <- markernames(brassica_genes_reduced_2, chr='A03')
chr3
chr4 <- markernames(brassica_genes_reduced_2, chr='A04')
chr4
chr5 <- markernames(brassica_genes_reduced_2, chr='A05')
chr5
chr6 <- markernames(brassica_genes_reduced_2, chr='A06')
chr6
chr7 <- markernames(brassica_genes_reduced_2, chr='A07')
chr7
chr8 <- markernames(brassica_genes_reduced_2, chr='A08')
chr8
chr9 <- markernames(brassica_genes_reduced_2, chr='A09')
chr9
# A10 is flipped, reverse marker order 
chr10 <- markernames(brassica_genes_reduced_2, chr='A10')
chr10 <- rev(chr10)
chr10

brassica_genes_reduced_2 <- switch.order(test, chr = "A10", order = chr10)
markernames(brassica_genes_reduced_2, chr='A10')

plotMap(brassica_genes_reduced_2, chr='A10') 
plot.rf(brassica_genes_reduced_2, chr='A10')


# save output progress
setwd("/Users/Cody_2/git.repos/brassica_genetic_map/Output")
write.cross(brassica_genes_reduced_2, format= "csvsr", 
	         filestem="snp_map_rqtl_Mbp_ref1.5_cross_output")


plotMap(brassica_genes_reduced_2) 
plot.rf(brassica_genes_reduced_2)

# Jan 5, 2015- Small Update to fix marker genetic position on A08.
brassica_fix <- read.cross("csvsr", genfile ="snp_map_rqtl_Mbp_ref1.5_cross_output_gen.csv", 
	                       phefile="snp_map_rqtl_Mbp_ref1.5_cross_output_phe.csv", 
	                       genotypes=c("AA","BB"), 
	                       na.strings=c("-"))

head(brassica_fix)
summary(brassica_fix)

# change to RIL population
class(brassica_fix)[1] <- "riself"
brassica_fix<- jittermap(brassica_fix)
brassica_fix

plotMap(brassica_fix) 
summaryMap(brassica_fix)
#         n.mar length ave.spacing max.spacing
# A01       163   94.0         0.6         3.1
# A02       119  121.3         1.0        36.6
# A03       207  135.6         0.7         4.6
# A04       100   63.3         0.6         2.1
# A05        84   91.1         1.1        26.3
# A06       158   93.9         0.6         2.9
# A07       160  106.4         0.7         3.1
# A08        77   80.7         1.1        29.1
# A09       216  154.7         0.7        14.5
# A10       167  104.2         0.6         2.1
# overall  1451 1045.1         0.7        36.6
plot.rf(brassica_fix)

# takes about a minute to estimate the map
?est.map
newmap <- est.map(brassica_fix, error.prob= 0.0005)
summaryMap(newmap)
#         n.mar length ave.spacing max.spacing
# A01       163   92.5         0.6         3.1
# A02       119  120.5         1.0        36.7
# A03       207  135.6         0.7         4.6
# A04       100   63.3         0.6         2.1
# A05        84   90.7         1.1        26.4
# A06       158   93.9         0.6         2.9
# A07       160  106.4         0.7         3.1
# A08        77   80.6         1.1        29.1
# A09       216  151.2         0.7        14.4
# A10       167  101.2         0.6         2.1
# overall  1451 1035.9         0.7        36.7
plotMap(newmap)

?shiftmap
# shift all starting positions of the map to 0
newmap <- shiftmap(newmap)

brassica_fix <- replace.map(brassica_fix, newmap)
plotMap(brassica_fix)
plot.rf(brassica_fix, chr = "A10")
pull.map(brassica_fix)

setwd("/Users/Cody_2/git.repos/brassica_genetic_map/Output")
write.cross(brassica_fix, format= "csvsr", 
	         filestem="snp_map_rqtl_Mbp_ref1.5.1_cross_output")


###########
# compare physical vs. genetic postions
###########
library(plyr)
library(ggplot2)
#reformat for plotting in ggplot
v2.1_map <- pull.map(brassica_fix, as.table=TRUE)
v2.1_map
v2.1_map$markernames <- rownames(v2.1_map)
rownames(v2.1_map) <- NULL
head(v2.1_map)

# split marker names
v2.1_map_description <- ldply(strsplit(as.character(v2.1_map$markernames),split="x"))
head(v2.1_map_description)
v2.1_map$genomic_pos <- v2.1_map_description[,2]

head(v2.1_map)
str(v2.1_map)
plot(v2.1_map$pos, v2.1_map$genomic_pos)

v2.1_map$genomic_pos <- as.numeric(v2.1_map$genomic_pos)
v2.1_map$genomic_pos <- (v2.1_map$genomic_pos)/1000000


A01_plot <- ggplot(v2.1_map, aes(x=pos, y=genomic_pos)) 
plot <- plot + geom_point(shape=1) + facet_grid(chr ~ . ) 
plot <- plot + xlab("Genetic Postition (cM)") + ylab("Physical Postion (Mbp)")
plot

#plot genomic and genetic postion of each marker against one another
setwd("/Users/Cody_2/git.repos/brassica_genetic_map/Output")
v2.1_A01 <- subset(v2.1_map, chr == "A01")
plot(v2.1_A01$pos, v2.1_A01$genomic_pos)

head(v2.1_A01)
str(v2.1_A01)
# A01
v2.1_A01$genomic_pos <- as.numeric(v2.1_A01$genomic_pos)
A01_plot <- ggplot(v2.1_A01, aes(x=pos, y=genomic_pos)) 
A01_plot <- A01_plot + geom_point(shape=1) + facet_grid(chr ~ . ) 
A01_plot <- A01_plot + xlab("Genetic Postition (cM)") + ylab("Physical Postion (Mbp)")
A01_plot
ggsave("A01_genetic_vs_physical_v2.1.pdf")

# A02
v2.1_A02 <- subset(v2.1_map, chr == "A02")
plot(v2.1_A02$pos, v2.1_A02$genomic_pos)

head(v2.1_A02)
str(v2.1_A02)

v2.1_A02$genomic_pos <- as.numeric(v2.1_A02$genomic_pos)
A02_plot <- ggplot(v2.1_A02, aes(x=pos, y=genomic_pos)) 
A02_plot <- A02_plot + geom_point(shape=1) + facet_grid(chr ~ . ) 
A02_plot <- A02_plot + xlab("Genetic Postition (cM)") + ylab("Physical Postion (Mbp)")
A02_plot
ggsave("A02_genetic_vs_physical_v2.1.pdf")

# A03
v2.1_A03 <- subset(v2.1_map, chr == "A03")
plot(v2.1_A03$pos, v2.1_A03$genomic_pos)

head(v2.1_A03)
str(v2.1_A03)

v2.1_A03$genomic_pos <- as.numeric(v2.1_A03$genomic_pos)
A03_plot <- ggplot(v2.1_A03, aes(x=pos, y=genomic_pos)) 
A03_plot <- A03_plot + geom_point(shape=1) + facet_grid(chr ~ . ) 
A03_plot <- A03_plot + xlab("Genetic Postition (cM)") + ylab("Physical Postion (Mbp)")
A03_plot
ggsave("A03_genetic_vs_physical_v2.1.pdf")

# A04
v2.1_A04 <- subset(v2.1_map, chr == "A04")
plot(v2.1_A04$pos, v2.1_A04$genomic_pos)

head(v2.1_A04)
str(v2.1_A04)

v2.1_A04$genomic_pos <- as.numeric(v2.1_A04$genomic_pos)
A04_plot <- ggplot(v2.1_A04, aes(x=pos, y=genomic_pos)) 
A04_plot <- A04_plot + geom_point(shape=1) + facet_grid(chr ~ . ) 
A04_plot <- A04_plot + xlab("Genetic Postition (cM)") + ylab("Physical Postion (Mbp)")
A04_plot
ggsave("A04_genetic_vs_physical_v2.1.pdf")

# A05
v2.1_A05 <- subset(v2.1_map, chr == "A05")
plot(v2.1_A05$pos, v2.1_A05$genomic_pos)

head(v2.1_A05)
str(v2.1_A05)

v2.1_A05$genomic_pos <- as.numeric(v2.1_A05$genomic_pos)
A05_plot <- ggplot(v2.1_A05, aes(x=pos, y=genomic_pos)) 
A05_plot <- A05_plot + geom_point(shape=1) + facet_grid(chr ~ . ) 
A05_plot <- A05_plot + xlab("Genetic Postition (cM)") + ylab("Physical Postion (Mbp)")
A05_plot
ggsave("A05_genetic_vs_physical_v2.1.pdf")

# A06
v2.1_A06 <- subset(v2.1_map, chr == "A06")
plot(v2.1_A06$pos, v2.1_A06$genomic_pos)

head(v2.1_A06)
str(v2.1_A06)

v2.1_A06$genomic_pos <- as.numeric(v2.1_A06$genomic_pos)
A06_plot <- ggplot(v2.1_A06, aes(x=pos, y=genomic_pos)) 
A06_plot <- A06_plot + geom_point(shape=1) + facet_grid(chr ~ . ) 
A06_plot <- A06_plot + xlab("Genetic Postition (cM)") + ylab("Physical Postion (Mbp)")
A06_plot
ggsave("A06_genetic_vs_physical_v2.1.pdf")

# A07
v2.1_A07 <- subset(v2.1_map, chr == "A07")
plot(v2.1_A07$pos, v2.1_A07$genomic_pos)

head(v2.1_A07)
str(v2.1_A07)

v2.1_A07$genomic_pos <- as.numeric(v2.1_A07$genomic_pos)
A07_plot <- ggplot(v2.1_A07, aes(x=pos, y=genomic_pos)) 
A07_plot <- A07_plot + geom_point(shape=1) + facet_grid(chr ~ . ) 
A07_plot <- A07_plot + xlab("Genetic Postition (cM)") + ylab("Physical Postion (Mbp)")
A07_plot
ggsave("A07_genetic_vs_physical_v2.1.pdf")

# A08
v2.1_A08 <- subset(v2.1_map, chr == "A08")
plot(v2.1_A08$pos, v2.1_A08$genomic_pos)

head(v2.1_A08)
str(v2.1_A08)

v2.1_A08$genomic_pos <- as.numeric(v2.1_A08$genomic_pos)
A08_plot <- ggplot(v2.1_A08, aes(x=pos, y=genomic_pos)) 
A08_plot <- A08_plot + geom_point(shape=1) + facet_grid(chr ~ . ) 
A08_plot <- A08_plot + xlab("Genetic Postition (cM)") + ylab("Physical Postion (Mbp)")
A08_plot
ggsave("A08_genetic_vs_physical_v2.1.pdf")

# A09
v2.1_A09 <- subset(v2.1_map, chr == "A09")
plot(v2.1_A09$pos, v2.1_A09$genomic_pos)

head(v2.1_A09)
str(v2.1_A09)

v2.1_A09$genomic_pos <- as.numeric(v2.1_A09$genomic_pos)
A09_plot <- ggplot(v2.1_A09, aes(x=pos, y=genomic_pos)) 
A09_plot <- A09_plot + geom_point(shape=1) + facet_grid(chr ~ . ) 
A09_plot <- A09_plot + xlab("Genetic Postition (cM)") + ylab("Physical Postion (Mbp)")
A09_plot
ggsave("A09_genetic_vs_physical_v2.1.pdf")

# A10
v2.1_A10 <- subset(v2.1_map, chr == "A10")
plot(v2.1_A10$pos, v2.1_A10$genomic_pos)

head(v2.1_A10)
str(v2.1_A10)

v2.1_A10$genomic_pos <- as.numeric(v2.1_A10$genomic_pos)
A10_plot <- ggplot(v2.1_A10, aes(x=pos, y=genomic_pos)) 
A10_plot <- A10_plot + geom_point(shape=1) + facet_grid(chr ~ . ) 
A10_plot <- A10_plot + xlab("Genetic Postition (cM)") + ylab("Physical Postion (Mbp)")
A10_plot
ggsave("A10_genetic_vs_physical_v2.1.pdf")

# end



















