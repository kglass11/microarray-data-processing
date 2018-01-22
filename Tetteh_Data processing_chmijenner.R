#####################################
###PROCESSING YOUR MICROARRAY data###
#####################################

###Install any packages you may need for this script
library(limma)
library(contrast)
library(beeswarm)
library(mixtools)
library(gplots)
library(ggplot2)
library(gcookbook)

###Generate specific data frames e.g. median minus background & mean minus background
#The column identifying the sample_id in slides_all.df id column 42
#The column identifying the foreground median in slides_all.df is 9
#The column identifying the background median in slides_all.df is 14
#The column identifying the median.df - the background in slides_all.df is 34 (here we ignore this column)

#Foreground
fore.df <- annotation_targets.df
for(i in 1:length(samples)){
  ite_out<-slides_all.df[which(slides_all.df[,42]==samples_unique[i]),9]
  fore.df<-cbind(fore.df,ite_out)
  colnames(fore.df)[length(colnames(fore.df))]<-samples_unique[i]
}

#Background
back.df <- annotation_targets.df
for(i in 1:length(samples)){
  ite_out<-slides_all.df[which(slides_all.df[,42]==samples_unique[i]),14]
  back.df<-cbind(back.df,ite_out)
  colnames(back.df)[length(colnames(back.df))]<-samples_unique[i]
}

###Generate matrices of the same data (that is, the same data with only a single data type)

fore.matrix <- as.matrix(fore.df[,6:ncol(fore.df)])
back.matrix <- as.matrix(back.df[,6:ncol(back.df)])

#Generate a corrected matrix, with background substracted from foreground
cor.matrix <- fore.matrix-back.matrix

#Adjust the corrected dataframe to get rid of zeros and normalise the values close to zero, using the limma pachage
######################


#Export this for reference
write.csv(cor.matrix, file="cor.matrix.csv")

###Assign target names to groups of your array targets to identify their 'type'

targets_blank = c(grep("BLANK", annotation_targets.df$Name))
targets_buffer = c(grep("buffer", annotation_targets.df$Name))
targets_ref = c(grep("REF", annotation_targets.df$Name))
targets_std = c(grep("Std", annotation_targets.df$Name))

###Begin quality control
#Check variability of control targets
#For this, the matrix must be transposed. We do this within the command using 't'

##1. Background
#To identify which slides, pads, and samples are significantly deviated, we need to calculate the mean and generate an arbitrary cut off 
#The cut off can be used to flag samples, and tells us whether deviation is universal, or specific to slides, pads, or samples

#Mean median background MFI for all data points
back_mean <- mean(back.matrix)
#SD median background MFI fro all data points
back_sd <- sd(back.matrix)
#Cut-off for deviation from mean
back_cutoff <- back_mean+(3*back_sd)

#Identify deviant samples/slides/pads
#
#Generate a vector of mean background intensity for every target for EACH sample (background person magnitude)
back_sample_mean <- colMeans(back.matrix)
#Vectors to identify normal and deviant samples
back_normal <- which(back_sample_mean<=back_cutoff)
back_deviant <- which(back_sample_mean>back_cutoff)
#Generate a table showing which slides, pads, samples have deviant background values
back_sample_deviant <- samples.df[back_deviant,]
write.csv(back_sample_deviant, file = "deviant_sample_background.csv")
#Coefficient of variation for all background datapoints, and all exlcuding deviant samples
back_cov_all <- sd(back.matrix)/mean(back.matrix)
back_cov_normal <- sd(back.matrix[, back_normal])/mean(back.matrix[, back_normal])

#Plots by slide/pad/sample
png(filename = paste("background_mfi_samples.tif"), width = 5.5, height = 10, units = "in", res = 600)
par(mfrow=c(3,1), mar = c(2, 3, 2.25, 0.5), oma = c(11.5, 0, 1, 0), bty = "o", 
    mgp = c(2, 0.5, 0), cex.main = 1.5, cex.axis = 0.6, cex.lab = 0.9, xpd=NA, las=1)
boxplot(t(back.matrix) ~ samples.df$slide_no, outcex=0.5,
        ylab="Background MFI", xlab="Slide", add=FALSE)
abline(h = back_cutoff, col = "red", lty = 2, lwd = 0.7, xpd=FALSE)
title(main = "Median background MFI by SLIDE/SUBARRAY/SAMPLE\n", adj=0)
boxplot(t(back.matrix) ~ samples.df$block_rep_1, outcex=0.5,
        ylab="Background MFI", xlab="Subarrays (1/2, 3/4, etc.)", add=FALSE, las=1)
abline(h = back_cutoff, col = "red", lty = 2, lwd = 0.7, xpd=FALSE)
boxplot(t(back.matrix) ~ samples.df$sample_id_unique, outcex=0.5,
        ylab="Background MFI", xlab="Samples", add=FALSE, las=2, cex.axis = 0.4, yaxt="n")
axis(2, cex.axis=0.6)
abline(h = back_cutoff, col = "red", lty = 2, lwd = 0.7, xpd=FALSE)

mtext(c(paste("Mean overall background MFI:", round(back_mean, digits=3))), side=1, cex=0.8, line=2, outer=TRUE, xpd=NA, adj=0)
mtext(c(paste("SD overall background MFI:", round(back_sd, digits=3))), side=1, cex=0.8, line=3.5, outer=TRUE, xpd=NA, adj=0)
mtext(c(paste("Cut-off overall background MFI:", round(back_cutoff, digits=3))), side=1, cex=0.8, line=5, outer=TRUE, xpd=NA, adj=0)
mtext(c(paste("CoV overall background MFI:", round(back_cov_all, digits=3))), side=1, cex=0.8, line=6.5, outer=TRUE, xpd=NA, adj=0)
mtext(c(paste("CoV overall background MFI (excl. deviant samples):", round(back_cov_normal, digits=3))), side=1, cex=0.8, line=8, outer=TRUE, xpd=NA, adj=0)

graphics.off()

#Identify deviant targets across all pads
#
#Mean background MFI for every sample for EACH target (background protein magnitude)
back_target_mean <- rowMeans(back.matrix)
#Mean background target magnitude for block 1, arranged in the order they are printed
back_target_mean_b1 = matrix(back_target_mean[annotation_targets.df$Block==1], nrow=max(annotation_targets.df$Row), ncol=max(annotation_targets.df$Column))
#Mean background target magnitude for block 2, arranged in the order they are printed
back_target_mean_b2 = matrix(back_target_mean[annotation_targets.df$Block==2], nrow=max(annotation_targets.df$Row), ncol=max(annotation_targets.df$Column))
#Are any background values for specific targets universally deviant, accross all pads?
back_target_deviant <- annotation_targets.df[back_target_mean>back_cutoff,]
back_target_deviant
write.csv(back_target_deviant, file = "deviant_target_background.csv")

#Plot background by target
png(filename = paste("background_mfi_targets.tif"), width = 15, height = 10, units = "in", res = 600)
par(mfrow=c(2,1), mar = c(7, 3, 2.25, 0.5), oma = c(6, 0, 1, 0), bty = "o", 
    mgp = c(2, 0.5, 0), cex.main = 1.5, cex.axis = 0.3, cex.lab = 0.9, xpd=NA, las=2)
boxplot(t(back.matrix), outcex=0.5, yaxt="n", ylab="Background MFI")
title(main="Background median MFI by ARRAY TARGET: all samples", adj=0)
abline(h = back_cutoff, col = "red", lty = 2, lwd = 0.7, xpd=FALSE)
axis(2, cex.axis=0.6)
mtext(c(paste("Targets 1-", index_target)), side=1, cex=0.9, xpd=NA, line=4.5, las=1)
boxplot(t(back.matrix[,back_normal]), outcex=0.5, yaxt="n", ylab="Background MFI")
title(main="Background median MFI by ARRAY TARGET: excl. deviant samples", adj=0)
axis(2, cex.axis=0.6)
mtext(c(paste("Targets 1-", index_target)), side=1, cex=0.9, xpd=NA, line=4.5, las=1)

mtext(c(paste("Mean overall background MFI:", round(back_mean, digits=3))), side=1, cex=0.8, line=0, outer=TRUE, xpd=NA, adj=0, las=1)
mtext(c(paste("SD overall background MFI:", round(back_sd, digits=3))), side=1, cex=0.8, line=1, outer=TRUE, xpd=NA, adj=0, las=1)
mtext(c(paste("Cut-off overall background MFI:", round(back_cutoff, digits=3))), side=1, cex=0.8, line=2, outer=TRUE, xpd=NA, adj=0, las=1)
mtext(c(paste("CoV overall background MFI:", round(back_cov_all, digits=3))), side=1, cex=0.8, line=3, outer=TRUE, xpd=NA, adj=0, las=1)
mtext(c(paste("CoV overall background MFI (excl. deviant samples):", round(back_cov_normal, digits=3))), side=1, cex=0.8, line=4, outer=TRUE, xpd=NA, adj=0, las=1)

graphics.off()

##2. Buffer
#Though some slides/pads may have high background, background correction may adjust this and make the data useable.
#The only way of identifying if this is true is by looking at the control spots by slide and pad.
#If the same samples have high controls, the background correction was not enough.
#Correction against plate buffers may be required...

#To identify which slides, pads, and samples are significantly deviated, we need to calculate the mean and generate an arbitrary cut off 
#The cut off can be used to flag samples, and tells us whether deviation is universal, or specific to slides, pads, or samples

#Mean median corrected MFI for all data points
cor_buffer_mean <- mean(cor.matrix[targets_buffer,])
#SD median corrected MFI for all data points
cor_buffer_sd <- sd(cor.matrix[targets_buffer,])
#Cut-off for deviation from mean
cor_cutoff <- cor_buffer_mean+(3*cor_buffer_sd)

#Identify deviant samples/slides/pads
#
#Generate a vector of mean buffer intensity for EACH sample (corrected person magnitude)
cor_buffer_sample_mean <- colMeans(cor.matrix[targets_buffer,])
#Vectors to identify normal and deviant samples
cor_normal <- which(cor_buffer_sample_mean<=cor_cutoff)
cor_deviant <- which(cor_buffer_sample_mean>cor_cutoff)
#Generate a table showing which slides, pads, samples have deviant corrected buffer values
cor_sample_deviant <- samples.df[cor_deviant,]
write.csv(cor_sample_deviant, file = "deviant_sample_corrected.csv")
#Coefficient of variation for all background datapoints, and all exlcuding deviant samples
cor_buffer_cov_all <- sd(cor.matrix[targets_buffer,])/mean(cor.matrix[targets_buffer,])
cor_buffer_cov_all_normal <- sd(cor.matrix[targets_buffer, cor_normal])/mean(cor.matrix[targets_buffer, cor_normal])

#Plots by slide/pad/sample
png(filename = paste("cor_mfi_all.tif"), width = 5.5, height = 10, units = "in", res = 600)
par(mfrow=c(3,1), mar = c(2, 3, 2.25, 0.5), oma = c(11.5, 0, 1, 0), bty = "o", 
    mgp = c(2, 0.5, 0), cex.main = 1.5, cex.axis = 0.6, cex.lab = 0.9, xpd=NA, las=1)
boxplot(t(cor.matrix[targets_buffer,]) ~ samples.df$slide_no, outcex=0.5,
        ylab="Corrected MFI", xlab="Slide", add=FALSE)
abline(h = cor_cutoff, col = "red", lty = 2, lwd = 0.7, xpd=FALSE)
title(main = "Median corrected MFI (buffer spots only)\n", adj=0)
boxplot(t(cor.matrix[targets_buffer,]) ~ samples.df$block_rep_1, outcex=0.5,
        ylab="Corrected MFI", xlab="Subarrays (1/2, 3/4, etc.)", add=FALSE, las=1)
abline(h = cor_cutoff, col = "red", lty = 2, lwd = 0.7, xpd=FALSE)
boxplot(t(cor.matrix[targets_buffer,]) ~ samples.df$sample_id_unique, outcex=0.5,
        ylab="Corrected MFI", xlab="Sample", add=FALSE, las=2, cex.axis = 0.4, yaxt="n")
axis(2, cex.axis=0.6)
abline(h = cor_cutoff, col = "red", lty = 2, lwd = 0.7, xpd=FALSE)

mtext(c(paste("Mean overall corrected buffer MFI:", round(cor_buffer_mean, digits=3))), side=1, cex=0.8, line=2, outer=TRUE, xpd=NA, adj=0)
mtext(c(paste("SD overall corrected buffer MFI:", round(cor_buffer_sd, digits=3))), side=1, cex=0.8, line=3.5, outer=TRUE, xpd=NA, adj=0)
mtext(c(paste("Cut-off overall corrected buffer MFI:", round(cor_cutoff, digits=3))), side=1, cex=0.8, line=5, outer=TRUE, xpd=NA, adj=0)
mtext(c(paste("CoV overall corrected buffer MFI:", round(cor_buffer_cov_all, digits=3))), side=1, cex=0.8, line=6.5, outer=TRUE, xpd=NA, adj=0)
mtext(c(paste("CoV overall corrected buffer MFI (excl. deviant samples):", round(cor_buffer_cov_all_normal, digits=3))), side=1, cex=0.8, line=8, outer=TRUE, xpd=NA, adj=0)

graphics.off()

#Identify deviant buffer targets across all pads
#
#Mean corrected MFI for every sample for EACH target (background protein magnitude)
cor_target_mean <- rowMeans(cor.matrix)
#Are any corrected buffer targets universally deviant, accross all pads?
temp <- which(cor_target_mean>cor_cutoff)
temp2 <- Reduce(intersect, list(targets_buffer, temp))
cor_buffer_deviant <- annotation_targets.df[temp2,]
cor_buffer_deviant
remove(temp, temp2)
write.csv(back_target_deviant, file = "deviant_buffer_corrected.csv")

#All slides
png(filename = paste("buffer_spots.tif"), width = 5, height = 4, units = "in", res = 600)
par(mar = c(5, 3, 2.25, 0.5), oma = c(0, 0, 0, 0), bty = "o", 
    mgp = c(2, 0.5, 0), cex.main = 1, cex.axis = 0.6, cex.lab = 0.7, xpd=NA, las=2)
boxplot(t(cor.matrix[targets_buffer,]),
        cex=0.5,
        ylab="Raw MFI")
abline(h = cor_cutoff, col = "red", lty = 2, lwd = 0.7, xpd=FALSE)
graphics.off()

#By slide
#mfrow is ordered by the number of rows, and columns of plots you want - so must be edited based on the number of slides
png(filename = paste("buffer_spots_slide.tif"), width = 5, height = 4, units = "in", res = 600)
par(mfrow=c(2,3), mar = c(4, 3, 1, 0.5), oma = c(1, 1, 1, 1), bty = "o", 
    mgp = c(2, 0.5, 0), cex.main = 1, cex.axis = 0.5, cex.lab = 0.7, xpd=NA, las=2)
for (i in 1:index_slide){
  boxplot(t(cor.matrix[targets_buffer,samples.df$slide_no==i]),
          ylab="Raw MFI",
          ylim=c(0,2000),
          add=FALSE, 
          cex=0.5,
          xpd=NA,
          main=c(paste("Slide", i)),
          adj=0)
  abline(h = cor_cutoff, col = "red", lty = 2, lwd = 0.7, xpd=FALSE)
}
graphics.off()

#Plot CoV by sample
cor_buffer_cov_sample <- c()

png(filename = paste("cov_buffer.tif"), width = 5, height = 4, units = "in", res = 600)
par(mar = c(4, 3, 1, 0.5), oma = c(1, 1, 1, 1), bty = "o", 
    mgp = c(2, 0.5, 0), cex.main = 1, cex.axis = 0.5, cex.lab = 0.7, xpd=NA, las=1)

for(i in 1:ncol(cor.matrix))
{
  temp <- sd(cor.matrix[targets_buffer,i])/mean(cor.matrix[targets_buffer,i])
  cor_buffer_cov_sample <- c(cor_buffer_cov_sample, temp)
}
plot(cor_buffer_cov_sample, ylab="CoV corrected buffer MFI", xlab="Sample", ylim=c(0,max(cor_buffer_cov_sample)))
remove(temp)
graphics.off()

#Check the number of values below 0, and convert them to 1
cor2.matrix <- cor.matrix
length(cor2.matrix)
sum(cor2.matrix<=0)
sum(cor2.matrix>0)
cat("Before baseline adjustment", sum(cor2.matrix<=0),"out of", length(cor2.matrix), "total targets have MFIs below zero in the corrected matrix")

#Simple method of baseline adjustment
for (i in 1:ncol(cor2.matrix)) 
{
  cor2.matrix[which(cor2.matrix[, i] < 1), i] <- 1
}
cat("After baseline adjustment", sum(cor2.matrix<=0),"out of", length(cor2.matrix), "total targets have MFIs below zero in the corrected matrix")

#Checking standards
plot(cor.matrix[targets_std,1])
axis(side=1, at=c(1:12), labels=rownames(cor.matrix[targets_std, ]))








######Stopped here in Ghana analysis
#So this is in wrong order, should be 
#1. Correction for zeros
#2. Log transform
#3. Correction for background (generating ratios)
#4. Normalisation against curves?

#Here, weve skipped this, so it is 
#1. Corrected for background
#2. Corrected for zeros

#For these figures, we will just remove the buffer from the corrected dataframes, rather than the log transformed ones
buffer_mean <- c()

for(i in 1:ncol(median.matrix))
{
  temp <- mean(median.matrix[buffer.probes,i])
  buffer_mean <- c(buffer_mean, temp)
}

#substract buffer mean from each sample to generate new intensity matrices
#remove highly variable samples from analyses
cor.median.matrix <- median.matrix
for(i in 1:ncol(cor.median.matrix))
{
  cor.median.matrix[,i] <- cor.median.matrix[,i]-buffer_mean[i]
}

#This is old - should really be corrected using Limma models
#Check the number of values below 0, and convert them to 1
#cor2.median.matrix <- cor.median.matrix
#for (i in 1:ncol(cor2.median.matrix)) 
#{
#  cor.median.matrix[which(cor2.median.matrix[, i] < 1), i] <- 1
#}
#Thus
cor.median.matrix <- as.data.frame(backgroundCorrect(as.matrix(cor.median.matrix), method = "normexp", offset = 50, normexp.method = "mle"))

cat("Numer of negative values in UNCORRECTED matrix", sum(median.matrix < 1), "out of", sum(cor2.median.matrix>0), "(", sum(median.matrix < 1)/sum(cor2.median.matrix>0)*100, "%)") 
cat("Numer of negative values in UNCORRECTED matrix", sum(cor.median.matrix < 1), "out of", sum(cor2.median.matrix>0), "(", sum(cor.median.matrix < 1)/sum(cor2.median.matrix>0)*100, "%)") 
cat("Numer of negative values in CORRECTED matrix", sum(cor.median.matrix < 1))

#Create unique rownames for proteins for combined matrices
ann.spot = median.df[1:(nrow(median.df)/2), 1:4]
ann.spot.unique <- rownames(ann.spot)
rownames(ann.spot) <- c(paste(ann.spot.unique, ann.spot$Name, sep = "_"))

#Make smaller corrected data frames
median.matrix.b1 <- median.matrix[which(median.df[,1]==1),]
median.matrix.b2 <- median.matrix[which(median.df[,1]==2),]
median.matrix.average <- (median.matrix.b1+median.matrix.b2)/2

rownames(median.matrix.b1) <- rownames(ann.spot)
rownames(median.matrix.b2) <- rownames(ann.spot)
rownames(median.matrix.average) <- rownames(ann.spot)

#Importing metadata
metadata.df = read.csv("All_samples.csv")
rownames(metadata.df) <- metadata.df$Sample
#metadata2.df <- metadata.df[samples_chr,]
metadata.df <- metadata.df[metadata.df$Sample %in% samples_chr, ]

population <- droplevels(population)

# Match colnames and pheno.df$sampleid
index = match(colnames(median.matrix.average), metadata.df$Sample)
metadata.df = metadata.df[index, ]

# Figures
age_breaks = c(0, 5, 15)
age_labels = c("0-5", "5-15")
age_2groups <- cut(metadata.df$Age, breaks= age_breaks, right= FALSE, labels= age_labels)
table(age_2groups)

png(filename = paste("buffer_cov.tif"), width = 3, height = 4, units = "in", res = 600)
par(mar = c(2.25, 2.25, 2.25, 0.5), oma = c(0, 0, 0, 0), bty = "o", 
    mgp = c(1.25, 0.5, 0), cex.main = 1, cex.axis = 0.6, cex.lab = 0.7, xpd=NA)

plot(cov_buffer, 
     ylab="CoV of buffer spots (n=8 per sample)", 
     xlab="Sample index (n=192)",
     col=ifelse(cov_buffer>2.5|cov_buffer<(-2.5), "red", "blue"))
legend(-15,-35,rmsamp_cov_chr, text.col ="red", adj = 0, pch=NA, cex = 0.6, bty="n")
graphics.off()

correl <- c()
for (i in 1:ncol(median.matrix.b1))
{
  correl_coef_temp <- cor.test(median.matrix.b1[,i], median.matrix.b2[,i])
  correl <- as.numeric(c(correl, correl_coef_temp[4]))
}  

png(filename = paste("block_correl.tif"), width = 3, height = 4, units = "in", res = 600)
par(mar = c(2.25, 2.25, 2.25, 0.5), oma = c(0, 0, 0, 0), bty = "o", 
    mgp = c(1.25, 0.5, 0), cex.main = 1, cex.axis = 0.6, cex.lab = 0.7, xpd=NA)

plot(correl, 
     ylab="Correlation coef. between spot repeats (n=180 per sample)", 
     xlab="Sample index (n=192)")
     #col=ifelse(cov_buffer>2.5|cov_buffer<(-2.5), "red", "blue"))
#legend(-15,-35,rmsamp_cov_chr, text.col ="red", adj = 0, pch=NA, cex = 0.6, bty="n")
graphics.off()

#Markers
png(filename = paste("markers_raw.tif"), width = 5, height = 3, units = "in", res = 600)
par(mfrow = c(1,3), mar = c(2.25, 2.25, 2.25, 0.5), oma = c(0, 0, 0, 0), bty = "o", 
    mgp = c(1.25, 0.5, 0), cex.main = 1, cex.axis = 0.6, cex.lab = 0.7, xpd=NA)

beeswarm(median.matrix.average[149,]~age_2groups, col="grey", cex=0.5, xlab="", ylab="")
boxplot(median.matrix.average[149,]~age_2groups, add=TRUE, col="transparent", title(main="AMA-1", adj=0), xlab="Age", ylab="Corrected log2 MFI")
beeswarm(median.matrix.average[147,]~age_2groups, col="grey", cex=0.5, xlab="", ylab="")
boxplot(median.matrix.average[147,]~age_2groups, add=TRUE, col="transparent", title(main="GLURP-R2", adj=0), xlab="Age", ylab="Corrected log2 MFI")
beeswarm(median.matrix.average[6,]~age_2groups, col="grey", cex=0.5, xlab="", ylab="")
boxplot(median.matrix.average[6,]~age_2groups, add=TRUE, col="transparent", title(main="ETRAMP-2(1)", adj=0), xlab="Age", ylab="Corrected log2 MFI")

graphics.off()

beeswarm(median.matrix.average[149,]~metadata.df$Population)

#Controls
blank.probes.average = c(grep("BLANK", ann.spot$Name))
buffer.probes.average = c(grep("buffer", ann.spot$Name))
ref.probes.average = c(grep("REF", ann.spot$Name))
std.probes.average = c(grep("Std", ann.spot$Name))

#loop doesnt work
for(i in 1:ncol(median.matrix.average))
{
plot(median.matrix.average[std.probes.average,1], xaxt='n', ylab="Raw MFI", ylim=c(0,30000))
  par(new=TRUE)
}
axis(side=1, at=c(1:6), labels=c("200ug/ml", "100ug/ml", "50ug/ml", "25ug/ml", "12.5ug/ml", "6.25ug/ml"), cex=0.2)

#Individual plots
png(filename = paste("markers_raw.tif"), width = 7, height = 3, units = "in", res = 600)
par(mfrow = c(1,5), mar = c(2.25, 2.25, 2.25, 0.5), oma = c(0, 0, 0, 0), bty = "o", 
    mgp = c(1.25, 0.5, 0), cex.main = 1, cex.axis = 0.5, cex.lab = 0.6, xpd=NA)

plot(median.matrix.average[std.probes.average,4], xaxt='n', ylab="Raw MFI", xlab="Standard curve (IgG[ug/ml])", ylim=c(0,30000))
axis(side=1, at=c(1:6), labels=c("200", "100", "50", "25", "12.5", "6.25"), cex=0.2, title(main = c(samples_chr[4]), adj=0))
plot(median.matrix.average[std.probes.average,5], xaxt='n', ylab="Raw MFI", xlab="Standard curve (IgG[ug/ml])", ylim=c(0,30000))
axis(side=1, at=c(1:6), labels=c("200", "100", "50", "25", "12.5", "6.25"), cex=0.2, title(main = c(samples_chr[5]), adj=0))
plot(median.matrix.average[std.probes.average,6], xaxt='n', ylab="Raw MFI", xlab="Standard curve (IgG[ug/ml])", ylim=c(0,30000))
axis(side=1, at=c(1:6), labels=c("200", "100", "50", "25", "12.5", "6.25"), cex=0.2, title(main = c(samples_chr[6]), adj=0))
plot(median.matrix.average[std.probes.average,7], xaxt='n', ylab="Raw MFI", xlab="Standard curve (IgG[ug/ml])", ylim=c(0,30000))
axis(side=1, at=c(1:6), labels=c("200", "100", "50", "25", "12.5", "6.25"), cex=0.2, title(main = c(samples_chr[7]), adj=0))
plot(median.matrix.average[std.probes.average,8], xaxt='n', ylab="Raw MFI", xlab="Standard curve (IgG[ug/ml])", ylim=c(0,30000))
axis(side=1, at=c(1:6), labels=c("200", "100", "50", "25", "12.5", "6.25"), cex=0.2, title(main = c(samples_chr[8]), adj=0))

graphics.off()

max(median.matrix.average[std.probes.average, ])
min(median.matrix.average[std.probes.average, ])

##########################################
#Convert values <0
#Cor3 here will be the data proceeding in the correct way
cor3.median.matrix <- median.matrix
for (i in 1:ncol(cor3.median.matrix)) 
{
  cor3.median.matrix[which(cor3.median.matrix[, i] < 1), i] <- 1
}

#Transform the database to log2
log2.median.matrix <- log2(cor3.median.matrix)

#Generate vector of mean buffer intensity for all samples for data correction
buffer_mean <- c()

for(i in 1:ncol(log2.median.matrix))
{
  temp <- mean(log2.median.matrix[buffer.probes,i])
  buffer_mean <- c(buffer_mean, temp)
}

#substract buffer mean from each sample to generate new intensity matrices
#remove highly variable samples from analyses
cor.log2.median.matrix <- log2.median.matrix
for(i in 1:ncol(cor.log2.median.matrix))
{
  cor.log2.median.matrix[,i] <- cor.log2.median.matrix[,i]-buffer_mean[i]
}

#Smaller dataframes
cor.log2.median.matrix.b1 <- cor.log2.median.matrix[which(median.df[,1]==1),]
cor.log2.median.matrix.b2 <- cor.log2.median.matrix[which(median.df[,1]==2),]
cor.log2.median.matrix.average <- (cor.log2.median.matrix.b1+cor.log2.median.matrix.b2)/2

rownames(cor.log2.median.matrix.b1) <- rownames(ann.spot)
rownames(cor.log2.median.matrix.b2) <- rownames(ann.spot)
rownames(cor.log2.median.matrix.average) <- rownames(ann.spot)

png(filename = paste("markers_log2_age.tif"), width = 5, height = 3, units = "in", res = 600)
par(mfrow = c(1,3), mar = c(2.25, 2.25, 2.25, 0.5), oma = c(0, 0, 0, 0), bty = "o", 
    mgp = c(1.25, 0.5, 0), cex.main = 1, cex.axis = 0.6, cex.lab = 0.7, xpd=NA)

beeswarm(cor.log2.median.matrix.average[149,]~age_2groups, col="grey", cex=0.5, xlab="", ylab="")
boxplot(cor.log2.median.matrix.average[149,]~age_2groups, add=TRUE, col="transparent", title(main="AMA-1", adj=0), xlab="Age", ylab="Corrected log2 MFI")
beeswarm(cor.log2.median.matrix.average[147,]~age_2groups, col="grey", cex=0.5, xlab="", ylab="")
boxplot(cor.log2.median.matrix.average[147,]~age_2groups, add=TRUE, col="transparent", title(main="GLURP-R2", adj=0), xlab="Age", ylab="Corrected log2 MFI")
beeswarm(cor.log2.median.matrix.average[6,]~age_2groups, col="grey", cex=0.5, xlab="", ylab="")
boxplot(cor.log2.median.matrix.average[6,]~age_2groups, add=TRUE, col="transparent", title(main="ETRAMP-2(1)", adj=0), xlab="Age", ylab="Corrected log2 MFI")

graphics.off()

png(filename = paste("markers_log2_country.tif"), width = 5, height = 3, units = "in", res = 600)
par(mfrow = c(1,3), mar = c(2.25, 2.25, 2.25, 0.5), oma = c(0, 0, 0, 0), bty = "o", 
    mgp = c(1.25, 0.5, 0), cex.main = 1, cex.axis = 0.6, cex.lab = 0.7, xpd=NA)

beeswarm(cor.log2.median.matrix.average[149,]~population, col="grey", cex=0.5, xlab="", ylab="")
boxplot(cor.log2.median.matrix.average[149,]~population, add=TRUE, col="transparent", title(main="AMA-1", adj=0), xlab="", ylab="Corrected log2 MFI")
beeswarm(cor.log2.median.matrix.average[147,]~population, col="grey", cex=0.5, xlab="", ylab="")
boxplot(cor.log2.median.matrix.average[147,]~population, add=TRUE, col="transparent", title(main="GLURP-R2", adj=0), xlab="", ylab="Corrected log2 MFI")
beeswarm(cor.log2.median.matrix.average[6,]~population, col="grey", cex=0.5, xlab="", ylab="")
boxplot(cor.log2.median.matrix.average[6,]~population, add=TRUE, col="transparent", title(main="ETRAMP-2(1)", adj=0), xlab="", ylab="Corrected log2 MFI")

graphics.off()

#heatmap
heat.data = log2.median.matrix[order(index, decreasing = TRUE), order(population, metadata.df$Age)]  
heat.data = as.matrix(heat.data)

png(filename = paste("heatmap.tif"),
    width = 5.5, height = 6, units = "in", res = 600)
heatmap(heat.data, Rowv = NA, Colv = NA, cexRow = 0.15, cexCol = 0.15)
graphics.off()

# Build a targets and nodna and purified proteins dataframes
##NOT DONE
blank.df <- log2.df[blank.probes, ]
target.df <- log2.df[target.probes, ]
ann.spot.target <- ann.spot[target.probes, ]
ivtt.df <- log2.df[nodna.probes, ]
ann.spot.ivtt <- ann.spot[nodna.probes, ]
pure.df <- cor.df[pure.probes, ]
ann.spot.pure <- ann.spot[pure.probes, ]

# Leaving out the purified proteins, we here generate the mean no dna levels for each subarray (within each sample)

##UNNECESSARY
#ivtt.mean1 = colMeans(ivtt.df[which(ann.spot.ivtt$"Array.Column" == 1), ])
#ivtt.mean2 = colMeans(ivtt.df[which(ann.spot.ivtt$"Array.Column" == 2), ])
#ivtt.mean3 = colMeans(ivtt.df[which(ann.spot.ivtt$"Array.Column" == 3), ])
#ivtt.mean4 = colMeans(ivtt.df[which(ann.spot.ivtt$"Array.Column" == 4), ])

# Normalize the data by subtracting the noDNA background by subarray
#norm.target.raw.df = target.df
#for (i in 1:4) {
#  for (j in 1:ncol(norm.target.raw.df)) {
#    norm.target.raw.df[which(ann.spot.target$"Array.Column" == i), j] = 
#      norm.target.raw.df[which(ann.spot.target$"Array.Column" == i), j] - 
#      get(paste("ivtt.mean", i, sep = ""))[j]
#  }
#}

#norm.ivtt.raw.df = ivtt.df
#for (i in 1:4) {
#  for (j in 1:ncol(norm.ivtt.raw.df)) {
#    norm.ivtt.raw.df[which(ann.spot.ivtt$"Array.Column" == i), j] = 
#      norm.ivtt.raw.df[which(ann.spot.ivtt$"Array.Column" == i), j] - 
#      get(paste("ivtt.mean", i, sep = ""))[j]
#  }
#}

#Removing values below 1/4 of the noDNA (very minimal numbers) which could affect mixture model based analysis
##UNNECESSARY because there has been no NoDNA subtraction, so all values still above 0

#norm.target.df <- norm.target.raw.df
#norm.ivtt.df <- norm.ivtt.raw.df
#norm.target.df[norm.target.df < -2] <- -2
#norm.ivtt.df[norm.ivtt.df < -2] <- -2

# Calculate the background threshold
##START HERE

back = as.matrix(blank.df)
background = mean(back) + 2 * sd(back)
seropos = median.df > background
colnames(seropos) = colnames(median.df)
reactivity = rowSums(seropos) > ncol(median.df) * 0.05
cat("There are", paste(sum(reactivity == TRUE)), "reative antigens out of",
    nrow(median.df), "printed antigens")
#
seropos.numeric <- seropos+0

# Calculate breadth
breadth = colSums(seropos)

#DETERMINING SEROPOSITIVITY
#First transpose the target matrix.

norm.target.transposed.df = as.data.frame(t(norm.target.df))
#If need to remove certain samples before continuing. 
#Here i'm just replacing one column with good data, then i'll remove the data afterwards. 
#This is just so the model runs. Only necessary with Cameroon only so far (353).
#norm.target.transposed.df [,353] = norm.target.transposed.df [,352]
#And - Gambia (222)
#norm.target.transposed.df [,222] = norm.target.transposed.df [,221]
#Burkina faso ()
#norm.target.transposed.df [,353] = norm.target.transposed.df [,352]

###Continue here
cutoff_list <- c()
restarts_list <-c()

for(i in 1:528){
  
  temp <-  normalmixEM(norm.target.transposed.df[,i], lambda = NULL, mu = NULL, sigma = NULL, maxrestarts=100000)
  if ( temp$mu[1] < temp$mu[2]) cutoff_temp <- temp$mu[1]+3*temp$sigma[1]
  if (temp$mu[1] > temp$mu[2])  cutoff_temp <- temp$mu[2]+3*temp$sigma[2]
  restarts_temp <- temp$restarts
  cutoff_list <- c(cutoff_list,cutoff_temp)
  restarts_list <- c(restarts_list,restarts_temp)

}

seropos_mix <- norm.target.transposed.df
for(i in 1:528)
{
  if(restarts_list[i]<=100000){
    
    seropos_mix[,i] <- as.numeric((norm.target.transposed.df[,i] > cutoff_list[i]))
    
  } else {
    
    seropos_mix[,i]<-rep(NA,dim(seropos_mix)[1])
    
  }  
}

###Now to remove samples with >=3 restarts
#rmSamp_restarts = restarts_list >=3
#seropos_mix2 = seropos_mix [, -(rmSamp_restarts$FALSE+1)]
#restarts_list <- as.data.frame(restarts_list)
#for(i in 1:528)
#seropos_mix[,193] = NA
#seropos_mix[,220] = NA
#seropos_mix[,222] = NA
#seropos_mix[,251] = NA
#seropos_mix[,262] = NA
#seropos_mix[,316] = NA
#seropos_mix[,447] = NA
#seropos_mix[,484] = NA
# And the one i duplicated
#seropos_mix[,353] = NA

#write.csv(seropos_mix, file="seropos_mix.csv")

#seropos_mix.transposed = t(seropos_mix)

###SEROPOSITIVITY for PURE probes only

pure.transposed.df = as.data.frame(t(pure.df))
#If need to remove certain samples before continuing. 
#norm.target.transposed.df [,353] = norm.target.transposed.df [,352]

###Continue here
pure_cutoff_list <- c()
pure_restarts_list <-c()

for(i in 1:22){
  
  temp <-  normalmixEM(pure.transposed.df[,i], lambda = NULL, mu = NULL, sigma = NULL, maxrestarts=100000)
  if ( temp$mu[1] < temp$mu[2]) cutoff_temp <- temp$mu[1]+3*temp$sigma[1]
  if (temp$mu[1] > temp$mu[2])  cutoff_temp <- temp$mu[2]+3*temp$sigma[2]
  restarts_temp <- temp$restarts
  pure_cutoff_list <- c(pure_cutoff_list,cutoff_temp)
  pure_restarts_list <- c(pure_restarts_list,restarts_temp)
  
}

pure_seropos_mix <- pure.transposed.df
for(i in 1:22)
{
  if(pure_restarts_list[i]<=100000){
    
    pure_seropos_mix[,i] <- as.numeric((pure.transposed.df[,i] > pure_cutoff_list[i]))
    
  } else {
    
    pure_seropos_mix[,i]<-rep(NA,dim(pure_seropos_mix)[1])
    
  }  
}

#Bad script section

#Alternative for assignsing sample IDs in slides_all.df based on order of samples among blocks specifically within attached excel file. Blocks must be ordered sequentially for this to work.

m <- 1

while (m <= sample_index) {
  
  for (i in 1:slide_index) {
    
    for (j in seq(from=1,to=block_index-1,by=2)) {
      
      slides_all.df$Sample[slides_all.df$slide_no==i&(slides_all.df$Block==j|slides_all.df$Block==j+1)] <- as.character(samples[m])
      
      m <- m + 1
      
    }
  }
}

