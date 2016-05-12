#########Merge Maps to fill in marker gaps################
#Last updated 2014_01_22
#by Cody Markelz, markelz@gmail.com
#F5 map created by Fede
#Cody map created with RNA-seq SNP markers and RQTL
#Marc map created with same set or RNA-seq SNP markers and a combination
###of RQTL and Onemap
################

setwd("/Users/Cody_2/git.repos/brassica_genetic_map/input")
library(qtl)

genotypes <- read.csv("B.rapa_RIL_GTs_fede_11_07.csv",na.strings=c("","-"))
head(genotypes)
summary(genotypes)
#get fede map into same format as RQTL
colnames(genotypes) <- sub("X","RIL_",colnames(genotypes))
head(genotypes)

dim(genotypes)
str(genotypes)

genotypes$chr <- as.factor(genotypes$chr)
genotypes$cm <- as.factor(genotypes$cm)

cols <- 4:163
genotypes[,cols][genotypes[,cols] == "0"] <- "AA" #R500
genotypes[,cols][genotypes[,cols] == "1"] <- "-" #het
genotypes[,cols][genotypes[,cols] == "2"] <- "BB" #IMB211
head(genotypes)

?gsub()
genotypes$chr <- gsub("([1-9])", "A0\\1", genotypes$chr)
genotypes$chr <- gsub("A010", "A10", genotypes$chr)
head(genotypes)
tail(genotypes)

names(genotypes)[1] <- paste("id")
names(genotypes)[2] <- paste("chr")
names(genotypes)[3] <- paste("cm")
head(genotypes)
dim(genotypes)


#infile cody map
#brassica_genes_reduced_4_gen.csv
cody_geno <- read.csv("brassica_genes_reduced_4_gen.csv")
head(cody_geno)

names(cody_geno)[2] <- paste("chr")
names(cody_geno)[3] <- paste("cm")
head(cody_geno)
dim(cody_geno)


cody_genotypes_merged <- merge(cody_geno, genotypes, all = T)
head(cody_genotypes_merged)
tail(cody_genotypes_merged)

dim(cody_genotypes_merged)
head(cody_genotypes_merged[,128:165])

cody_genotypes_merged <- cody_genotypes_merged[,1:127]
head(cody_genotypes_merged)
tail(cody_genotypes_merged)
dim(cody_genotypes_merged)


names(cody_genotypes_merged)[2] <- ""
names(cody_genotypes_merged)[3] <- ""
head(cody_genotypes_merged)
tail(cody_genotypes_merged)
write.table(cody_genotypes_merged, file= "cody_old_merge_rqtl.csv", row.names = FALSE, col.names = TRUE,
	           sep = ",")

#cody and fede map merged
#markers behave poorly
library(qtl)

merged_map <- read.cross("csvsr", genfile ="cody_old_merge_rqtl.csv", 
	                       phefile="phenotype.csv", 
	                       genotypes=c("AA","BB"), 
	                       na.strings=c("NA","-"))

class(merged_map)[1] <- "riself"
merged_map <- jittermap(merged_map)
merged_map <- est.rf(merged_map)
merged_map <- drop.markers(merged_map, "FLC1aE")

print(dup <- findDupMarkers(merged_map, exact.only=FALSE))
head(dup)
merged_map <- drop.markers(merged_map, unlist(dup))
totmar(merged_map)
merged_map_2 <- merged_map
plotMap(merged_map)



mn2_merged_map <- markernames(merged_map, chr="A02")
mn2_merged_map

mn3 <- markernames(merged_map, chr="A03")
mn4 <- markernames(merged_map, chr="A04")
mn5 <- markernames(merged_map, chr="A05")
mn6 <- markernames(merged_map, chr="A06")
mn7 <- markernames(merged_map, chr="A07")
mn8 <- markernames(merged_map, chr="A08")
mn9 <- markernames(merged_map, chr="A09")
mn10 <- markernames(merged_map, chr="A10")
mn2
plotMap(merged_map, chr = "A02", show.marker.names=TRUE)

merged_map_3 <- merged_map_2
merged_map_4 <- merged_map_2
#drop markers on chromosome 2
merged_map_3 <- drop.markers(merged_map_3, 
	                     c("BRMS046","pX123bH", "fito338a", "pW130bX",
	                       "pX139eH",  "pW208aH", "pX110cX", "fito378a",
	                       "pW235aX", "pW250aH", "fito473",  "EZ3bH",   
	                       "pW135aH",  "pX128bX", "pW205bH", "pW176aH",  
	                       "pX107aX", "pW213aE", "pW249aX",  "pX124bX")) 
plotMap(merged_map_3, chr = "A02", show.marker.names=TRUE)

merged_map_4 <- drop.markers(merged_map_4, 
	                     c("fito529", "FLC2aE", "BRMS001",  "fito118b",
	                       "pW251aE" ,     "BRMS026a", "BRMS044", "fito338b",     
	                       "pW227cE","pW250aH", "fito473",  "EZ3bH",   
	                       "pW135aH",  "pX128bX", "pW205bH", "pW176aH",  
	                       "pX107aX", "pW213aE", "pW249aX",  "pX124bX",
	                       "pX123bH", "pW235aX" )) 
plotMap(merged_map_4, chr = "A02", show.marker.names=TRUE)
merged_map_4 <- orderMarkers(merged_map_4, chr=c('A02'), 
	                        window=5, use.ripple=TRUE, maxit=4000, 
	                        error.prob=0.1)
merged_map_4_new <- est.map(merged_map_4, chr = "A02", error.prob= 0.001, omit.noninformative=TRUE)
plotMap(merged_map_4)



newmap_3 <- est.map(merged_map_3, chr = "A02", error.prob= 0.001, omit.noninformative=TRUE)
summaryMap(newmap_3)

plotMap(newmap_3, merged_map_3, chr = "A02")
plotMap(newmap_3)
summaryMap(merged_map_3)

merged_map_3 <- orderMarkers(merged_map_3, chr=c('A02'), 
	                        window=7, use.ripple=TRUE, maxit=4000, 
	                        error.prob=0.1)
merged_map_3 <- est.rf(merged_map_3)
plot.rf(merged_map_3, chr='A02')
plotMap(merged_map_3, chr="A02")

?pull.rf
rf <- pull.rf(merged_map_4, chr = "A02")
lod <- pull.rf(merged_map_4, chr = "A02", what = "lod")
mn2 <- markernames(merged_map_4, chr="A02")
mn2

summaryMap(merged_map_4)
mn2_x <- markernames(merged_map, chr="A02")
mn2_x
#AO2_2670143
plot(rf, mn2[1], bandcol="gray70", ylim=c(0,1), alternate.chrid=TRUE)
plot(rf, mn2[2], bandcol="gray70", ylim=c(0,1), alternate.chrid=TRUE)
plot(rf, mn2[3], bandcol="gray70", ylim=c(0,1), alternate.chrid=TRUE)
plot(rf, mn2[4], bandcol="gray70", ylim=c(0,1), alternate.chrid=TRUE)
plot(rf, mn2[5], bandcol="gray70", ylim=c(0,1), alternate.chrid=TRUE)

plot(rf, mn2[107], bandcol="gray70", ylim=c(0,1), alternate.chrid=TRUE)
plot(rf, mn2[16], bandcol="gray70", ylim=c(0,1), alternate.chrid=TRUE)






#split
plot(rf, mn2[6], bandcol="gray70", ylim=c(0,1), alternate.chrid=TRUE)
plot(rf, mn2[32], bandcol="gray70", ylim=c(0,1), alternate.chrid=TRUE)
plot(rf, mn2[35], bandcol="gray70", ylim=c(0,1), alternate.chrid=TRUE)
plot(rf, mn2[40], bandcol="gray70", ylim=c(0,1), alternate.chrid=TRUE)
plot(rf, mn2[31], bandcol="gray70", ylim=c(0,1), alternate.chrid=TRUE)
#to marker 40, 
plot(rf, mn2[29], bandcol="gray70", ylim=c(0,1), alternate.chrid=TRUE)

plot(rf, mn2[100], bandcol="gray70", ylim=c(0,1), alternate.chrid=TRUE)
plotMap(merged_map, merged_map_3, chr = "A02")


?est.map
?est.rf
newmap <- est.map(merged_map, chr = "A02", error.prob= 0.1, omit.noninformative=TRUE)
summaryMap(newmap)
summaryMap(merged_map)
?plotMap
merged_map <- replace.map(merged_map, newmap)

mn2 <- markernames(merged_map, chr="A02")
mn2
plotMap(merged_map, chr= "A02", show.marker.names=TRUE)
merged_map <- drop.markers(merged_map, c("fito356","pW249dX","fito389","fito367a", ))
newmap <- est.map(merged_map, chr = "A02", error.prob= 0.1, omit.noninformative=TRUE)
merged_map <- replace.map(merged_map, newmap)

summaryMap(merged_map)
summaryMap(newmap)

plotMap(newmap)

#play with marker orders
#this is as if I am starting over
merged_map <- replace.map(merged_map,newmap)

?ripple
?orderMarkers()
merged_map <- orderMarkers(merged_map, chr=c('A02'), 
	                        window=2, use.ripple=TRUE, maxit=4000, 
	                        error.prob=0.0001)
plot.rf(merged_map, chr='A02')


merged_map_2 <- orderMarkers(merged_map_2, chr = "A01",
	                        window=2, use.ripple=TRUE, maxit=4000, 
	                        error.prob=0.0001)

rip_1 <- ripple(merged_map_2, chr = "A01", window=2, method = "likelihood")
summary(rip_1)
merged_map_2 <- switch.order(merged_map_2, chr = "A01", rip_1[2,])
plotMap(merged_map_2, chr = "A01")


plotMap(merged_map_2, chr = "A02")
merged_map_2 <- orderMarkers(merged_map_2, chr = "A02",
	                        window=2, use.ripple=TRUE, maxit=4000, 
	                        error.prob=0.0001)

rip_2 <- ripple(merged_map_2, chr = "A02", window=2, method = "likelihood")
summary(rip_2)
merged_map_2 <- switch.order(merged_map_2, chr = "A02", rip_2[2,])
plotMap(merged_map_2, chr = "A02")


plotMap(merged_map_2, chr = "A08")
merged_map_2 <- orderMarkers(merged_map_2, chr = "A08",
	                        window=2, use.ripple=TRUE, maxit=4000, 
	                        error.prob=0.0001)

rip_8 <- ripple(merged_map_2, chr = "A08", window=2, method = "likelihood")
summary(rip_8)
merged_map_2 <- switch.order(merged_map_2, chr = "A08", rip_8[2,])
?plotMap
plotMap(merged_map_2, chr = "A08", show.marker.names= TRUE)




####################
####################
####################
####cody marc map merged analysis
####################
####################
####################
#infile cody map
#brassica_genes_reduced_4_gen.csv
cody_geno <- read.csv("brassica_genes_reduced_4_gen.csv")
head(cody_geno)

names(cody_geno)[2] <- paste("chr")
names(cody_geno)[3] <- paste("cm")
head(cody_geno)
dim(cody_geno)

#infile marc map
marc_geno <- read.csv("snp_map_rqtl_MARC.csv")
head(marc_geno)


names(marc_geno)[2] <- paste("chr")
names(marc_geno)[3] <- paste("cm")
head(marc_geno)
dim(marc_geno)

cody_marc_merged <- merge(cody_geno, marc_geno, all = TRUE)
head(cody_marc_merged)
dim(cody_marc_merged)
tail(cody_marc_merged)
cody_marc_merged[1:20,]

#remove duplicated markers
cody_marc_merged_2 <- cody_marc_merged[which( !duplicated(cody_marc_merged[,c("id")]) ), ]
head(cody_marc_merged_2)
dim(cody_marc_merged_2)
tail(cody_marc_merged_2)


names(cody_marc_merged_2)[2] <- ""
names(cody_marc_merged_2)[3] <- ""
head(cody_marc_merged_2)
tail(cody_marc_merged_2)
write.table(cody_marc_merged_2, file= "cody_marc_merge_rqtl.csv", row.names = FALSE, col.names = TRUE,
	           sep = ",")
library(qtl)

merged_map <- read.cross("csvsr", genfile ="cody_marc_merge_rqtl.csv", 
	                       phefile="phenotype.csv", 
	                       genotypes=c("AA","BB"), 
	                       na.strings=c("NA","-", "AB"))

class(merged_map)[1] <- "riself"
merged_map <- jittermap(merged_map)
merged_map <- est.rf(merged_map)
summary(merged_map)
summaryMap(merged_map)
plotMap(merged_map, chr= "A08", show.marker.names = TRUE)
newmap <- est.map(merged_map, error.prob= 0.1)
plotMap(merged_map, newmap)


