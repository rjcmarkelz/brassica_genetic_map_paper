# Following Karl Broman's Guide Here: LINK



###
### Create genetic map
###

library(qtl)
setwd("/Users/Cody_2/git.repos/brassica_genetic_map_paper/input")

b_map <- read.cross("csvsr", genfile ="snp_map_rqtl_Mbp_ref1.5.csv", 
	                       phefile="brock_2010_pheno.csv", 
	                       genotypes=c("AA","BB"), 
	                       na.strings=c("NA","AB"))

head(b_map)
summary(b_map)

# change to RIL population
class(b_map)[1] <- "riself"
b_map <- jittermap(b_map)
b_map

# takes about a minute to estimate the map
newmap <- est.map(b_map, error.prob= 0.0005)
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

plotMap(newmap) 
plot.rf(b_map)  

#find potential duplicates and drop them
dup <- findDupMarkers(b_map, exact.only = FALSE)
head(dup)
length(dup)
# 57

totmar(b_map)
# 1527

b_map <- drop.markers(b_map, unlist(dup))
totmar(b_map)
# 1465

#identify markers with possible segregation distortion
gt <- geno.table(b_map)
seg_dist <- gt[gt$P.value < 0.05/totmar(b_map),]
seg_dist
dim(seg_dist)

# 33 all on chromosome 3, will not drop these, but will make note in paper.
#              chr missing AA  BB      P.value
# A03x19849963 A03       0 38  86 1.628658e-05
# A03x19928814 A03       0 36  88 3.015843e-06
# A03x20027810 A03       0 37  87 7.117887e-06
# A03x20139956 A03       0 36  88 3.015843e-06 ...... truncated output

# reestimate map with dropped duplicate markers
newmap <- est.map(b_map, error.prob = 0.0005)
summaryMap(newmap)
plotMap(newmap) # saved _2
b_map <- replace.map(b_map, newmap)
plot.rf(b_map) # saved _2



# now remove markers that obviously do not belong
# make a copy so I can go back and not have to reestimate the map
b_map_red <- b_map
rffull <- pull.rf(b_map_red)

#############
# warm up with A08
#############
plot.rf(b_map_red, chr = 'A08')
rf8 <- pull.rf(b_map_red, chr = 'A08')
chr8 <- markernames(b_map_red, chr = 'A08')
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
b_map_red <- drop.markers(b_map_red, c("A08x4875175",  "A08x4952134"))


b_map_red <- orderMarkers(b_map_red, chr = c('A08'), 
	                        window = 5, use.ripple = TRUE, maxit = 4000, 
	                        error.prob = 0.0001)
plotMap(b_map_red, chr ='A08') 
plot.rf(b_map_red, chr = 'A08')
rf8 <- pull.rf(b_map_red, chr = 'A08')
chr8 <- markernames(b_map_red, chr = 'A08')
chr8

#still something funny
plot(rffull, chr8[1], bandcol = "gray70", ylim = c(0,1), alternate.chrid = TRUE)
plot(rffull, chr8[2], bandcol = "gray70", ylim = c(0,1), alternate.chrid = TRUE)
plot(rffull, chr8[3], bandcol = "gray70", ylim = c(0,1), alternate.chrid = TRUE)
# markers "A08x1515279"  "A08x1485026" may belong on chromosome 2

b_map_red <- drop.markers(b_map_red, c("A08x1515279",  "A08x1485026"))

setseed(1234)
b_map_red <- orderMarkers(b_map_red, chr = c('A08'), 
	                        window = 5, use.ripple = TRUE, maxit = 4000, 
	                        error.prob = 0.0001)
plotMap(b_map_red, chr = 'A08') 
plot.rf(b_map_red, chr = 'A08')


# take a look at total map
plotMap(b_map_red)

#############
# Onto A05
#############
plot.rf(b_map_red, chr='A05')
rf5 <- pull.rf(b_map_red, chr = 'A05')
chr5 <- markernames(b_map_red, chr='A05')
chr5

plot(rf5, chr5[22], bandcol = "gray70", ylim = c(0,1), alternate.chrid = TRUE)
plot(rf5, chr5[23], bandcol = "gray70", ylim = c(0,1), alternate.chrid = TRUE)
plot(rf5, chr5[24], bandcol = "gray70", ylim = c(0,1), alternate.chrid = TRUE)
plot(rf5, chr5[25], bandcol = "gray70", ylim = c(0,1), alternate.chrid = TRUE)
plot(rf5, chr5[26], bandcol = "gray70", ylim = c(0,1), alternate.chrid = TRUE)
plot(rf5, chr5[27], bandcol = "gray70", ylim = c(0,1), alternate.chrid = TRUE)
plot(rf5, chr5[28], bandcol = "gray70", ylim = c(0,1), alternate.chrid = TRUE)

plot(rffull, chr5[22], bandcol = "gray70", ylim = c(0,1), alternate.chrid = TRUE)
plot(rffull, chr5[23], bandcol = "gray70", ylim = c(0,1), alternate.chrid = TRUE)

plot(rffull, chr5[24], bandcol = "gray70", ylim = c(0,1), alternate.chrid = TRUE)
plot(rffull, chr5[25], bandcol = "gray70", ylim = c(0,1), alternate.chrid = TRUE)
plot(rffull, chr5[26], bandcol = "gray70", ylim = c(0,1), alternate.chrid = TRUE)
plot(rffull, chr5[27], bandcol = "gray70", ylim = c(0,1), alternate.chrid = TRUE)

# this batch of bad markers belongs on chr A08
# "A05x8149908" "A05x8245839"  "A05x8277002"  "A05x8327143"
b_map_red <- drop.markers(b_map_red, c("A05x8149908", "A05x8245839",  "A05x8277002",  "A05x8327143"))

set.seed(1256)
b_map_red <- orderMarkers(b_map_red, chr = c('A05'), 
	                        window = 5, use.ripple = TRUE, maxit = 4000, 
	                        error.prob = 0.0001)
plotMap(b_map_red, chr = 'A05') 
plot.rf(b_map_red, chr = 'A05')

# find the next target
plotMap(b_map_red)

#############
# A02
#############
plot.rf(b_map_red, chr = 'A02')
rf2 <- pull.rf(b_map_red, chr  =  'A02')
chr2 <- markernames(b_map_red, chr = 'A02')
chr2
# marker "A02x21152172" is strange and belongs on A06
plot(rf2, chr2[76], bandcol = "gray70", ylim = c(0,1), alternate.chrid = TRUE)
plot(rf2, chr2[77], bandcol = "gray70", ylim = c(0,1), alternate.chrid = TRUE)
plot(rf2, chr2[78], bandcol = "gray70", ylim = c(0,1), alternate.chrid = TRUE)
plot(rf2, chr2[79], bandcol = "gray70", ylim = c(0,1), alternate.chrid = TRUE)

plot(rffull, chr2[78], bandcol = "gray70", ylim = c(0,1), alternate.chrid = TRUE)

# marker "A02x23373611" is strange and belongs on chromosome A10
plot(rf2, chr2[91], bandcol = "gray70", ylim = c(0,1), alternate.chrid = TRUE)
plot(rf2, chr2[92], bandcol = "gray70", ylim = c(0,1), alternate.chrid = TRUE)
plot(rf2, chr2[93], bandcol = "gray70", ylim = c(0,1), alternate.chrid = TRUE)

plot(rffull, chr2[93], bandcol = "gray70", ylim = c(0,1), alternate.chrid = TRUE)

b_map_red <- drop.markers(b_map_red, c("A02x21152172", "A02x23373611"))

# may need to come back to 2, it is still really large ~200 cM
set.seed(1465)
b_map_red <- orderMarkers(b_map_red, chr = c('A02'), 
	                        window = 5, use.ripple = TRUE, maxit = 4000, 
	                        error.prob = 0.0001)
plotMap(b_map_red, chr = 'A02') 
plot.rf(b_map_red, chr = 'A02')

plot.rf(b_map_red, chr='A02')

rffull <- pull.rf(b_map_red)
rf2 <- pull.rf(b_map_red, chr = 'A02')
chr2 <- markernames(b_map_red, chr = 'A02')
chr2

plot(rf2, chr2[20], bandcol = "gray70", ylim = c(0,1), alternate.chrid = TRUE)
plot(rffull, chr2[70], bandcol = "gray70", ylim = c(0,1), alternate.chrid = TRUE)

# just needed to reorder the markers on the top of A02
set.seed(165)
b_map_red <- orderMarkers(b_map_red, chr = c('A02'), 
	                        window = 2,  maxit = 4000, 
	                        error.prob = 0.0001)
plotMap(b_map_red, chr = 'A02') 
plot.rf(b_map_red, chr = 'A02')


#############
# A01
#############
plot.rf(b_map_red, chr = 'A01')
rf1 <- pull.rf(b_map_red, chr  =  'A01')
chr1 <- markernames(b_map_red, chr = 'A01')
chr1

# marker "A01x8900522" belongs on A09
plot(rf1, chr1[58], bandcol = "gray70", ylim = c(0,1), alternate.chrid = TRUE)
plot(rf1, chr1[59], bandcol = "gray70", ylim = c(0,1), alternate.chrid = TRUE)
plot(rffull, chr1[58], bandcol = "gray70", ylim = c(0,1), alternate.chrid = TRUE)

# marker "A01x10301489" belongs on A07
plot(rf1, chr1[79], bandcol = "gray70", ylim = c(0,1), alternate.chrid = TRUE)
plot(rffull, chr1[79], bandcol = "gray70", ylim = c(0,1), alternate.chrid = TRUE)

# marker "A01x11450244" belongs on A04
plot(rf1, chr1[84], bandcol = "gray70", ylim = c(0,1), alternate.chrid = TRUE)
plot(rffull, chr1[84], bandcol = "gray70", ylim = c(0,1), alternate.chrid = TRUE)

# marker "A01x16479326" belongs on A02
plot(rf1, chr1[92], bandcol = "gray70", ylim = c(0,1), alternate.chrid = TRUE)
plot(rffull, chr1[92], bandcol = "gray70", ylim = c(0,1), alternate.chrid = TRUE)

b_map_red <- drop.markers(b_map_red, c("A01x8900522", "A01x10301489", "A01x11450244", "A01x16479326"))

# need to run twice with same seed to reproduce marker order
# ripple to save one state, and then ripple again to clean up
set.seed(163234)
?orderMarkers
b_map_red <- orderMarkers(b_map_red, chr = c('A01'), 
	                        window = 2, use.ripple = FALSE, maxit = 4000, 
	                        error.prob = 0.0001)
plotMap(b_map_red, chr = 'A01') 

plot.rf(b_map_red, chr = 'A01')

# save output progress
setwd("/Users/Cody_2/git.repos/brassica_genetic_map_paper/Output")
write.cross(b_map_red, format= "csvsr", filestem="snp_map_rqtl_Mbp_ref1.5_cross_output")

plotMap(b_map_red) 
plot.rf(b_map_red)


#############
# A09 ripple marker order
#############
plotMap(b_map_red, chr = 'A09') 
plot.rf(b_map_red, chr = 'A09')

set.seed(16374)
b_map_red <- orderMarkers(b_map_red, chr = c('A09'), 
	                        window = 2, use.ripple = TRUE, maxit = 4000, 
	                        error.prob = 0.0001)

plotMap(b_map_red, chr = 'A09') 
plot.rf(b_map_red, chr = 'A09')

set.seed(16374)
b_map_red <- orderMarkers(b_map_red, chr = c('A09'), 
	                        window = 2, use.ripple = TRUE, maxit = 4000, 
	                        error.prob = 0.0001)

plotMap(b_map_red, chr = 'A09') 
plot.rf(b_map_red, chr = 'A09')

#############
# A10 ripple marker order
#############
plotMap(b_map_red, chr='A10') 
plot.rf(b_map_red, chr = 'A10')

set.seed(163778)
b_map_red <- orderMarkers(b_map_red, chr = c('A10'), 
	                        window = 2, use.ripple = TRUE, maxit = 4000, 
	                        error.prob = 0.0001)

plotMap(b_map_red, chr='A10') 
plot.rf(b_map_red, chr = 'A10')

# save output progress
setwd("/Users/Cody_2/git.repos/brassica_genetic_map_paper/Output")
write.cross(b_map_red, format= "csvsr", filestem = "snp_map_rqtl_Mbp_ref1.5_cross_output")

##########################
##########################
# Fine tuning
##########################
##########################
summaryMap(b_map_red)
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
# Back to A08
#############
plot.rf(b_map_red, chr='A08')
rf8 <- pull.rf(b_map_red, chr = 'A08')
chr8 <- markernames(b_map_red, chr='A08')
chr8
set.seed(1234)
b_map_red <- orderMarkers(b_map_red, chr=c('A08'), 
	                        window=2, use.ripple = FALSE, maxit=4000, 
	                        error.prob=0.0001)
plotMap(b_map_red, chr='A08') 
plot.rf(b_map_red, chr='A08')


#############
# Back to A05
#############
rf5 <- pull.rf(b_map_red, chr = 'A05')
chr5 <- markernames(b_map_red, chr='A05')
chr5
plot.rf(b_map_red, chr='A05')
plotMap(b_map_red, chr='A05') 

set.seed(1234)
b_map_red <- orderMarkers(b_map_red, chr=c('A05'), 
	                        window=2, use.ripple = FALSE, maxit=4000, 
	                        error.prob=0.0001)
plotMap(b_map_red, chr='A05') 
plot.rf(b_map_red, chr='A05')


# final plots
plotMap(b_map_red) 
plot.rf(b_map_red)

# do the marker orders make sense?
chr1 <- markernames(b_map_red, chr='A01')
chr1 # reverse marker order

chr2 <- markernames(b_map_red, chr='A02')
chr2 # reverse marker order

chr3 <- markernames(b_map_red, chr='A03')
chr3 # okay

chr4 <- markernames(b_map_red, chr='A04')
chr4 # okay

chr5 <- markernames(b_map_red, chr='A05')
chr5 # okay

chr6 <- markernames(b_map_red, chr='A06')
chr6 # okay

chr7 <- markernames(b_map_red, chr='A07')
chr7 # okay

chr8 <- markernames(b_map_red, chr='A08')
chr8 # reverse marker order

chr9 <- markernames(b_map_red, chr='A09')
chr9 # reorder within, or potential inversion?
# 37 - 52

chr10 <- markernames(b_map_red, chr='A10')
chr10 # reverse marker order, potential inversion?

# reorder if necessary
chr1 <- markernames(b_map_red, chr='A01')
chr1 <- rev(chr1)
chr1

b_map_red <- switch.order(b_map_red, chr = "A01", order = chr1)
markernames(b_map_red, chr='A01')

plotMap(b_map_red, chr='A01') 
plot.rf(b_map_red, chr='A01')

chr2 <- markernames(b_map_red, chr='A02')
chr2 <- rev(chr2)
chr2

b_map_red <- switch.order(b_map_red, chr = "A02", order = chr2)
markernames(b_map_red, chr='A02')

plotMap(b_map_red, chr='A02') 
plot.rf(b_map_red, chr='A02')

chr8 <- markernames(b_map_red, chr='A08')
chr8 <- rev(chr8)
chr8

b_map_red <- switch.order(b_map_red, chr = "A08", order = chr8)
markernames(b_map_red, chr='A08')

plotMap(b_map_red, chr='A08') 
plot.rf(b_map_red, chr='A08')


chr9 <- markernames(b_map_red, chr='A09')
chr9 
chr9[37:52] <- rev(chr9[37:52])
chr9

b_map_red <- switch.order(b_map_red, chr = "A09", order = chr9)
markernames(b_map_red, chr='A09')

plotMap(b_map_red, chr='A09') 
plot.rf(b_map_red, chr='A09')
b_map_red <- est.rf(b_map_red)


chr10 <- markernames(b_map_red, chr='A10')
chr10 <- rev(chr10)
chr10

b_map_red <- switch.order(b_map_red, chr = "A10", order = chr10)
markernames(b_map_red, chr='A10')

plotMap(b_map_red, chr='A10') 
plot.rf(b_map_red, chr='A10')

# make sure all start at 0 cM
b_map_red <- shiftmap(b_map_red)

# save output
setwd("/Users/Cody_2/git.repos/brassica_genetic_map_paper/Output")
write.cross(b_map_red, format= "csvsr", filestem="snp_map_rqtl_Mbp_ref1.5_cross_output")

# final checks and save files
pdf('plotMap_new_output.pdf')
plotMap(b_map_red) 

pdf('plot.rf_new_output.pdf')
plot.rf(b_map_red) 
dev.off()

# end



















