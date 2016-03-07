###
### Create genetic map
###
### Last updated March 6, 2016
### With fixed marker order, 

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
#         n.mar length ave.spacing max.spacing
# A01       175  387.2         2.2       276.4
# A02       125 1236.6        10.0       966.8
# A03       216  135.6         0.6         4.6
# A04       107  110.1         1.0        46.8
# A05        88  270.0         3.1       178.5
# A06       169  166.4         1.0        37.5
# A07       165  125.1         0.8        11.0
# A08        81  262.8         3.3       137.7
# A09       224  325.8         1.5        49.6
# A10       172  121.4         0.7        21.2
# overall  1522 3140.9         2.1       966.8


plotMap(newmap) # saved as genetic_map_no_rearrangement
plot.rf(b_map)  # saved as recombination_heat_map_no_rearrangement

#find potential duplicates and drop them
dup <- findDupMarkers(b_map, exact.only = FALSE)
head(dup)
length(dup)
# 40

totmar(b_map)
# 1522

b_map <- drop.markers(b_map, unlist(dup))
totmar(b_map)
# 1482

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
# A01
#############
plot.rf(b_map_red, chr = 'A01')
plotMap(b_map_red, chr = 'A01')
rf1 <- pull.rf(b_map_red, chr  =  'A01')
chr1 <- markernames(b_map_red, chr = 'A01')
chr1

# ripple to save one state, and then ripple again to clean up
set.seed(163234)
# ripple four times, looks great
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
# A02
#############
plot.rf(b_map_red, chr = 'A02')
rf2 <- pull.rf(b_map_red, chr  =  'A02')
chr2 <- markernames(b_map_red, chr = 'A02')
chr2

set.seed(1465) 
# ripple a few times and everthing looks good
b_map_red <- orderMarkers(b_map_red, chr = c('A02'), 
	                        window = 5, use.ripple = TRUE, maxit = 4000, 
	                        error.prob = 0.0001)
plotMap(b_map_red, chr = 'A02') 
plot.rf(b_map_red, chr = 'A02')

#############
# A03
#############
plot.rf(b_map_red, chr = 'A03')
rf3 <- pull.rf(b_map_red, chr  =  'A03')
chr3 <- markernames(b_map_red, chr = 'A03')
chr3

set.seed(1465) 
# ripple a few times and everthing looks good
b_map_red <- orderMarkers(b_map_red, chr = c('A03'), 
	                        window = 5, use.ripple = TRUE, maxit = 4000, 
	                        error.prob = 0.0001)
plotMap(b_map_red, chr = 'A03') 
plot.rf(b_map_red, chr = 'A03')

#############
# A04
#############
plot.rf(b_map_red, chr = 'A04')
rf4 <- pull.rf(b_map_red, chr  =  'A04')
chr4 <- markernames(b_map_red, chr = 'A04')
chr4

# ripple a few times
b_map_red <- orderMarkers(b_map_red, chr = c('A04'), 
	                        window = 5, use.ripple = TRUE, maxit = 4000, 
	                        error.prob = 0.0001)
plotMap(b_map_red, chr = 'A04') 
plot.rf(b_map_red, chr = 'A04')

# save output progress
setwd("/Users/Cody_2/git.repos/brassica_genetic_map_paper/Output")
write.cross(b_map_red, format= "csvsr", filestem="snp_map_rqtl_Mbp_ref1.5_cross_output")

#############
# A05
#############
plot.rf(b_map_red, chr='A05')
rf5 <- pull.rf(b_map_red, chr = 'A05')
chr5 <- markernames(b_map_red, chr='A05')
chr5

# running ripple a 3 times gets to final output
b_map_red <- orderMarkers(b_map_red, chr = c('A05'), 
	                        window = 5, use.ripple = TRUE, maxit = 4000, 
	                        error.prob = 0.0001)
plotMap(b_map_red, chr = 'A05') 
plot.rf(b_map_red, chr = 'A05')

# find the next target
plotMap(b_map_red)

#############
# A06
#############
plot.rf(b_map_red, chr='A06')
rf6 <- pull.rf(b_map_red, chr = 'A06')
chr6 <- markernames(b_map_red, chr='A06')
chr6

b_map_red <- orderMarkers(b_map_red, chr = c('A06'), 
	                        window = 5, use.ripple = TRUE, maxit = 4000, 
	                        error.prob = 0.0001)
plotMap(b_map_red, chr = 'A06') 
plot.rf(b_map_red, chr = 'A06')

# save output progress
setwd("/Users/Cody_2/git.repos/brassica_genetic_map_paper/Output")
write.cross(b_map_red, format= "csvsr", filestem="snp_map_rqtl_Mbp_ref1.5_cross_output")

#############
# A07
#############
plot.rf(b_map_red, chr='A07')
rf7 <- pull.rf(b_map_red, chr = 'A07')
chr7 <- markernames(b_map_red, chr='A07')
chr7

b_map_red <- orderMarkers(b_map_red, chr = c('A07'), 
	                        window = 5, use.ripple = TRUE, maxit = 4000, 
	                        error.prob = 0.0001)
plotMap(b_map_red, chr = 'A07') 
plot.rf(b_map_red, chr = 'A07')

# save output progress
setwd("/Users/Cody_2/git.repos/brassica_genetic_map_paper/Output")
write.cross(b_map_red, format= "csvsr", filestem="snp_map_rqtl_Mbp_ref1.5_cross_output")

#############
# A08 
#############
plot.rf(b_map_red, chr = 'A08')
rf8 <- pull.rf(b_map_red, chr = 'A08')
chr8 <- markernames(b_map_red, chr = 'A08')
chr8

# a few ripples
b_map_red <- orderMarkers(b_map_red, chr = c('A08'), 
	                        window = 5, use.ripple = TRUE, maxit = 4000, 
	                        error.prob = 0.0001)
plotMap(b_map_red, chr ='A08') 
plot.rf(b_map_red, chr = 'A08')
rf8 <- pull.rf(b_map_red, chr = 'A08')
chr8 <- markernames(b_map_red, chr = 'A08')
chr8

# save output progress
setwd("/Users/Cody_2/git.repos/brassica_genetic_map_paper/Output")
write.cross(b_map_red, format= "csvsr", filestem="snp_map_rqtl_Mbp_ref1.5_cross_output")

#############
# A09 
#############
plotMap(b_map_red, chr = 'A09') 
plot.rf(b_map_red, chr = 'A09')

#ripple a few times
b_map_red <- orderMarkers(b_map_red, chr = c('A09'), 
	                        window = 2, use.ripple = TRUE, maxit = 4000, 
	                        error.prob = 0.0001)
plotMap(b_map_red, chr = 'A09') 
plot.rf(b_map_red, chr = 'A09')

#############
# A10 
#############
plotMap(b_map_red, chr='A10') 
plot.rf(b_map_red, chr = 'A10')

# ripple marker 
b_map_red <- orderMarkers(b_map_red, chr = c('A10'), 
	                        window = 2, use.ripple = TRUE, maxit = 4000, 
	                        error.prob = 0.0001)
plotMap(b_map_red, chr='A10') 
plot.rf(b_map_red, chr = 'A10')

# save output progress
setwd("/Users/Cody_2/git.repos/brassica_genetic_map_paper/Output")
write.cross(b_map_red, format= "csvsr", filestem = "snp_map_rqtl_Mbp_ref1.5_cross_output")

summaryMap(b_map_red)
#         n.mar length ave.spacing max.spacing
# A01       170  103.2         0.6         3.6
# A02       123  121.0         1.0        36.6
# A03       207  135.6         0.7         4.6
# A04       102   63.6         0.6         2.1
# A05        87   91.0         1.1        26.3
# A06       159   93.9         0.6         2.9
# A07       163  106.3         0.7         3.1
# A08        81   80.5         1.0        29.1
# A09       220  153.0         0.7        13.7
# A10       170  104.0         0.6         2.5
# overall  1482 1052.2         0.7        36.6


# final plots
plotMap(b_map_red) 
plot.rf(b_map_red)

# do the marker orders make sense?
chr1 <- markernames(b_map_red, chr='A01')
chr1 # okay

chr2 <- markernames(b_map_red, chr='A02')
chr2 # okay

chr3 <- markernames(b_map_red, chr='A03')
chr3 # okay

chr4 <- markernames(b_map_red, chr='A04')
chr4 # reverse marker order

chr5 <- markernames(b_map_red, chr='A05')
chr5 # okay

chr6 <- markernames(b_map_red, chr='A06')
chr6 # reverse marker order

chr7 <- markernames(b_map_red, chr='A07')
chr7 # reverse marker order

chr8 <- markernames(b_map_red, chr='A08')
chr8 # okay

chr9 <- markernames(b_map_red, chr='A09')
chr9 # reverse marker order

chr10 <- markernames(b_map_red, chr='A10')
chr10 # reverse marker order

############
# reverse A04, A06, A07, A09, A10
############

# A04
chr4 <- markernames(b_map_red, chr='A04')
chr4 <- rev(chr4)
chr4

b_map_red <- switch.order(b_map_red, chr = "A04", order = chr4)
markernames(b_map_red, chr='A04')

plotMap(b_map_red, chr='A04') 
plot.rf(b_map_red, chr='A04')

# A06
chr6 <- markernames(b_map_red, chr='A06')
chr6 <- rev(chr6)
chr6

b_map_red <- switch.order(b_map_red, chr = "A06", order = chr6)
markernames(b_map_red, chr='A06')

plotMap(b_map_red, chr='A06') 
plot.rf(b_map_red, chr='A06')

# A07
chr7 <- markernames(b_map_red, chr='A07')
chr7 <- rev(chr7)
chr7

b_map_red <- switch.order(b_map_red, chr = "A07", order = chr7)
markernames(b_map_red, chr='A07')

plotMap(b_map_red, chr='A07') 
plot.rf(b_map_red, chr='A07')

# A09
chr9 <- markernames(b_map_red, chr='A09')
chr9 <- rev(chr9)
chr9 

b_map_red <- switch.order(b_map_red, chr = "A09", order = chr9)
markernames(b_map_red, chr='A09')

plotMap(b_map_red, chr='A09') 
plot.rf(b_map_red, chr='A09')

# A10
chr10 <- markernames(b_map_red, chr='A10')
chr10 <- rev(chr10)
chr10

b_map_red <- switch.order(b_map_red, chr = "A10", order = chr10)
markernames(b_map_red, chr='A10')

plotMap(b_map_red, chr='A10') 
plot.rf(b_map_red, chr='A10')

# make sure all start at 0 cM
b_map_red <- shiftmap(b_map_red)
b_map_red
# save output
setwd("/Users/Cody_2/git.repos/brassica_genetic_map_paper/Output")
write.cross(b_map_red, format= "csvsr", filestem="snp_map_rqtl_Mbp_ref1.5_cross_output")

# final checks and save files
dev.new()
pdf('plotMap_new_output.pdf')
plotMap(b_map_red) 

pdf('plot.rf_new_output.pdf')
plot.rf(b_map_red) 
dev.off()

# end



















