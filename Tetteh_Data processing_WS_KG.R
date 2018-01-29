#####################################
###PROCESSING YOUR MICROARRAY data###
#####################################

###Install any packages you may need for this script. Go to Tools, install packages, or run
#install.packages(c("contrast", "beeswarm", "mixtools", "gplots", "ggplot2", "gcookbook"))

##Install limma
## try http:// if https:// URLs are not supported
#source("https://bioconductor.org/biocLite.R")
#biocLite("limma")

###Load packages needed for this script
library(limma)
library(contrast)
library(beeswarm)
library(mixtools)
library(gplots)
library(ggplot2)
library(gcookbook)

###SELECT RELEVANT DATA###

###Generate specific data frames e.g. median foreground and background
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

#Generate matrices of the same data (that is, the same data with only a single data type)

fore.matrix <- as.matrix(fore.df[,7:ncol(fore.df)])
back.matrix <- as.matrix(back.df[,7:ncol(back.df)])

###DATA CORRECTION###

### Perform background correction using limma function and normexp method.
cor.matrix <- backgroundCorrect.matrix(fore.matrix, back.matrix, method = "normexp", offset = 50, normexp.method = "mle")

#Export this for reference
write.csv(cor.matrix, file="corrected.matrix.csv")



###Assign target names to groups of your array targets to identify their 'type'

targets_blank = c(grep("BLANK", annotation_targets.df$Name))
targets_buffer = c(grep("buffer", annotation_targets.df$Name))
targets_ref = c(grep("REF", annotation_targets.df$Name))
targets_std = c(grep("Std", annotation_targets.df$Name))
targets_allcontrol = c(targets_blank, targets_buffer, targets_ref, targets_std)

###Assign sample type names, to idenfy control and test samples
# KG - this is creating empty vectors right now, the samples aren't previously typed
# to be test or control
samples_test <- samples.df$sample_type=="test"
samples_control <- samples.df$sample_type=="control"

###Remove bad spots from subsequent analysis

#First, we need to tell R which spots are the highest spots on the plate
#the spots printed immediatley after these can be made to have null MFIs
#This call will create a vector idetifying all ref spots, plus all std spots called "Std 1"
high_targets <- c(targets_ref, grep("Std 1", row.names(annotation_targets.df)))
#To identify the spots to disinclude, we need to identify the position of the next spot in the print run.
#Because block 2 is printed before block 1, we can subtract an entire block, minus 1 row (e.g. 180-12=168)
#THIS WILL CHANGE DEPENDING ON YOUR EXACT PLATE PLAN - SO ALTER THE FIGURE ACCORDINGLY. 
high_targets_disinclude <- ifelse(high_targets>=180, high_targets-168, high_targets+180)
high_targets_disinclude<-high_targets_disinclude[which(high_targets_disinclude<=336)] 

###Convert the spots to be disincluded to NAs in background corrected data
cor2.matrix <- cor.matrix
cor2.matrix[high_targets_disinclude, ] <- NA

###QUALITY CONTROL###

###Check variability of control targets
#For this, the matrix must be transposed. We do this within the command using 't'
# KG - don't need transposed matrices for background QC...

### 1. Background
#To identify which slides, pads, and samples are significantly deviated, we need to calculate the mean and generate an arbitrary cut off 
#The cut off can be used to flag samples, and tells us whether deviation is universal, or specific to slides, pads, or samples

#Mean of median background MFI for all data points
back_mean <- mean(back.matrix)
#SD median background MFI fro all data points
back_sd <- sd(back.matrix)
#Cut-off for deviation from mean
back_cutoff <- back_mean+(3*back_sd)

###Identify deviant samples/slides/pads

#Generate a vector of mean background intensity for every target for EACH sample (background person magnitude)
back_sample_mean <- colMeans(back.matrix)
#Vectors to identify normal and deviant samples
back_normal <- which(back_sample_mean<=back_cutoff)
back_deviant <- which(back_sample_mean>back_cutoff)
#Generate a table showing which slides, pads, samples have deviant background values
back_sample_deviant <- samples.df[back_deviant,]
back_sample_deviant
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

###Identify deviant targets across all pads

#Mean background MFI for every sample for EACH target (background protein magnitude)
back_target_mean <- rowMeans(back.matrix)
#Mean background target magnitude for block 1, arranged in the order they are printed
back_target_mean_b1 = matrix(back_target_mean[annotation_targets.df$Block==1], nrow=max(annotation_targets.df$Row), ncol=max(annotation_targets.df$Column))
#Mean background target magnitude for block 2, arranged in the order they are printed
back_target_mean_b2 = matrix(back_target_mean[annotation_targets.df$Block==2], nrow=max(annotation_targets.df$Row), ncol=max(annotation_targets.df$Column))
#Are any background values for specific targets universally deviant, accross all pads?
back_target_deviant <- annotation_targets.df[back_target_mean>back_cutoff,]
back_target_deviant

#Plot background by target
png(filename = paste("background_mfi_targets.tif"), width = 15, height = 10, units = "in", res = 600)
par(mfrow=c(2,1), mar = c(7, 3, 2.25, 0.5), oma = c(6, 0, 1, 0), bty = "o", 
    mgp = c(2, 0.5, 0), cex.main = 1.5, cex.axis = 0.3, cex.lab = 0.9, xpd=NA, las=2)
boxplot(t(back.matrix), outcex=0.5, yaxt="n", ylab="Background MFI")
title(main="Background MFI by ARRAY TARGET: all samples", adj=0)
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

#Cov sample (all samples and all targets)
png(filename = paste("cov_back_sample.tif"), width = 10, height = 4, units = "in", res = 600)
par(mfrow=c(1,2), mar = c(4, 3, 1, 0.5), oma = c(1, 1, 1, 1), bty = "o", 
    mgp = c(2, 0.5, 0), cex.main = 1, cex.axis = 0.5, cex.lab = 0.7, xpd=NA, las=1)

back_cov_sample <- c()
for(i in 1:ncol(back.matrix))
{
  temp <- sd(back.matrix[,i])/mean(back.matrix[,i])
  back_cov_sample <- c(back_cov_sample, temp)
}
plot(back_cov_sample, ylab="CoV background MFI", xlab="Sample", ylim=c(0,max(back_cov_sample)))
remove(temp)

back_cov_target <- c()
for(i in 1:nrow(back.matrix))
{
  temp <- sd(back.matrix[i,])/mean(back.matrix[i,])
  back_cov_target <- c(back_cov_target, temp)
}
plot(back_cov_target, ylab="CoV background MFI", xlab="Target", ylim=c(0,max(back_cov_target)))
remove(temp)
graphics.off()

### 2. Buffer
#Though some slides/pads may have high background, background correction may adjust this and make the data useable.
#The only way of identifying if this is true is by looking at the control spots by slide and pad.
#If the same samples have high controls, the background correction was not enough.
#Correction against plate buffers may be required...

#To identify which slides, pads, and samples are significantly deviated, we need to calculate the mean and generate an arbitrary cut off 
#The cut off can be used to flag samples, and tells us whether deviation is universal, or specific to slides, pads, or samples

#Mean median corrected MFI for all data points
cor_buffer_mean <- mean(cor2.matrix[targets_buffer,], na.rm = TRUE)
#SD median corrected MFI for all data points
cor_buffer_sd <- sd(cor2.matrix[targets_buffer,], na.rm = TRUE)
#Cut-off for deviation from mean
cor_cutoff <- cor_buffer_mean+(3*cor_buffer_sd)

###Identify deviant samples/slides/pads

#Generate a vector of mean buffer intensity for EACH sample (corrected person magnitude)
cor_buffer_sample_mean <- colMeans(cor2.matrix[targets_buffer,], na.rm = TRUE)
#Vectors to identify normal and deviant samples
cor_normal <- which(cor_buffer_sample_mean<=cor_cutoff)
cor_deviant <- which(cor_buffer_sample_mean>cor_cutoff)
#Generate a table showing which slides, pads, samples have deviant corrected buffer values
cor_sample_deviant <- samples.df[cor_deviant,]
cor_sample_deviant
write.csv(cor_sample_deviant, file = "deviant_sample_corrected.csv")
#Coefficient of variation for all background datapoints, and all exlcuding deviant samples
cor_buffer_cov_all <- cor_buffer_sd/cor_buffer_mean
cor_buffer_cov_all_normal <- sd(cor2.matrix[targets_buffer, cor_normal], na.rm = TRUE)/mean(cor2.matrix[targets_buffer, cor_normal], na.rm = TRUE)

#Plots by slide/pad/sample
png(filename = paste("cor_mfi_samples.tif"), width = 5.5, height = 10, units = "in", res = 600)
par(mfrow=c(3,1), mar = c(2, 3, 2.25, 0.5), oma = c(11.5, 0, 1, 0), bty = "o", 
    mgp = c(2, 0.5, 0), cex.main = 1.5, cex.axis = 0.6, cex.lab = 0.9, xpd=NA, las=1)
boxplot(t(cor2.matrix[targets_buffer,]) ~ samples.df$slide_no, outcex=0.5,
        ylab="Corrected MFI", xlab="Slide", add=FALSE)
abline(h = cor_cutoff, col = "red", lty = 2, lwd = 0.7, xpd=FALSE)
title(main = "Corrected MFI (buffer only) by SLIDE/SUBARRAY/SAMPLE\n", adj=0)
boxplot(t(cor2.matrix[targets_buffer,]) ~ samples.df$block_rep_1, outcex=0.5,
        ylab="Corrected MFI", xlab="Subarrays (1/2, 3/4, etc.)", add=FALSE, las=1)
abline(h = cor_cutoff, col = "red", lty = 2, lwd = 0.7, xpd=FALSE)
boxplot(t(cor2.matrix[targets_buffer,]) ~ samples.df$sample_id_unique, outcex=0.5,
        ylab="Corrected MFI", xlab="Sample", add=FALSE, las=2, cex.axis = 0.4, yaxt="n")
axis(2, cex.axis=0.6)
abline(h = cor_cutoff, col = "red", lty = 2, lwd = 0.7, xpd=FALSE)

mtext(c(paste("Mean overall corrected buffer MFI:", round(cor_buffer_mean, digits=3))), side=1, cex=0.8, line=2, outer=TRUE, xpd=NA, adj=0)
mtext(c(paste("SD overall corrected buffer MFI:", round(cor_buffer_sd, digits=3))), side=1, cex=0.8, line=3.5, outer=TRUE, xpd=NA, adj=0)
mtext(c(paste("Cut-off overall corrected buffer MFI:", round(cor_cutoff, digits=3))), side=1, cex=0.8, line=5, outer=TRUE, xpd=NA, adj=0)
mtext(c(paste("CoV overall corrected buffer MFI:", round(cor_buffer_cov_all, digits=3))), side=1, cex=0.8, line=6.5, outer=TRUE, xpd=NA, adj=0)
mtext(c(paste("CoV overall corrected buffer MFI (excl. deviant samples):", round(cor_buffer_cov_all_normal, digits=3))), side=1, cex=0.8, line=8, outer=TRUE, xpd=NA, adj=0)

graphics.off()

###Identify deviant buffer targets across all pads

#Mean corrected MFI for every sample for EACH target (background protein magnitude)
# !!!error - temp2 not found, leads to cor_buffer_deviant not found, fix this later if we want to include this
cor_target_mean <- rowMeans(cor2.matrix)
#Are any corrected buffer targets universally deviant, accross all pads?
temp <- Reduce(intersect, list(targets_buffer, which(cor_target_mean>cor_cutoff)))
cor_buffer_deviant <- annotation_targets.df[temp2,]
cor_buffer_deviant
remove(temp)

#All slides
png(filename = paste("buffer_spots.tif"), width = 5, height = 4, units = "in", res = 600)
par(mar = c(5, 3, 2.25, 0.5), oma = c(0, 0, 0, 0), bty = "o", 
    mgp = c(2, 0.5, 0), cex.main = 1, cex.axis = 0.6, cex.lab = 0.7, xpd=NA, las=2)
boxplot(t(cor2.matrix[targets_buffer,]),
        cex=0.5,
        ylab="Corrected MFI")
abline(h = cor_cutoff, col = "red", lty = 2, lwd = 0.7, xpd=FALSE)
graphics.off()

#By slide
#mfrow is ordered by the number of rows, and columns of plots you want - so must be edited based on the number of slides
png(filename = paste("buffer_spots_slide.tif"), width = 5, height = 4, units = "in", res = 600)
par(mfrow=c(2,3), mar = c(4, 3, 1, 0.5), oma = c(1, 1, 1, 1), bty = "o", 
    mgp = c(2, 0.5, 0), cex.main = 1, cex.axis = 0.5, cex.lab = 0.7, xpd=NA, las=2)
for (i in 1:index_slide){
  boxplot(t(cor2.matrix[targets_buffer,samples.df$slide_no==i]),
          ylab="Corrected MFI",
          ylim=c(0,2000),
          add=FALSE, 
          cex=0.5,
          xpd=NA,
          main=c(paste("Slide", i)),
          adj=0)
  abline(h = cor_cutoff, col = "red", lty = 2, lwd = 0.7, xpd=FALSE)
}
graphics.off()

#Plot CoV of buffer by sample
cor_buffer_cov_sample <- c()

png(filename = paste("cov_buffer.tif"), width = 5, height = 4, units = "in", res = 600)
par(mar = c(4, 3, 1, 0.5), oma = c(1, 1, 1, 1), bty = "o", 
    mgp = c(2, 0.5, 0), cex.main = 1, cex.axis = 0.5, cex.lab = 0.7, xpd=NA, las=1)

for(i in 1:ncol(cor2.matrix))
{
  temp <- sd(cor2.matrix[targets_buffer,i], na.rm = TRUE)/mean(cor2.matrix[targets_buffer,i], na.rm = TRUE)
  cor_buffer_cov_sample <- c(cor_buffer_cov_sample, temp)
}
plot(cor_buffer_cov_sample, ylab="CoV corrected buffer MFI", xlab="Sample", ylim=c(0,max(cor_buffer_cov_sample)))
remove(temp)
graphics.off()

### Mean values for each target ordered so it could be a heatmap - this is not currently log transformed or normalized data

#Check average corrected values for each target for all individuals
cor_target_mean <- rowMeans(cor2.matrix)
#Mean background target magnitude for block 1, arranged in the order they are printed
cor_target_mean_b1 = t(matrix(round(cor_target_mean[annotation_targets.df$Block==1], digits=2), nrow=max(annotation_targets.df$Column), ncol=max(annotation_targets.df$Row)))
#Mean background target magnitude for block 2, arranged in the order they are printed
cor_target_mean_b2 = t(matrix(round(cor_target_mean[annotation_targets.df$Block==2], digits=2), nrow=max(annotation_targets.df$Column), ncol=max(annotation_targets.df$Row)))
cor_target_mean_b1b2 <- rbind(cor_target_mean_b1, cor_target_mean_b2)
#Annotation plate maps
annotation_target_b1 <- t(matrix(annotation_targets.df$target_id_unique [annotation_targets.df$Block==1], nrow = max(annotation_targets.df$Column), ncol=max(annotation_targets.df$Row)))
annotation_target_b2 <- t(matrix(annotation_targets.df$target_id_unique [annotation_targets.df$Block==2], nrow = max(annotation_targets.df$Column), ncol=max(annotation_targets.df$Row)))
annotation_target_b1b2 <- rbind(annotation_target_b1, annotation_target_b2)
#Write csv, which can be presented as a heatmap
write.csv(cbind(cor_target_mean_b1b2, annotation_target_b1b2), file="corrected_target_mean.csv")
remove(cor_target_mean, cor_target_mean_b1, cor_target_mean_b2,cor_target_mean_b1b2, annotation_target_b1, annotation_target_b2, annotation_target_b1b2)

### Log transform the data (base 2)
log.cor.matrix <- log2(cor2.matrix)

### Normalization ###

# remove this part once checked whole script! make sure to replace later with earlier buffer values
###Create new buffer summary based on only good spots, background corrected, not log-transformed 
## checked, these values are now giving the same values as generated above in buffer section 
cor2_buffer_mean <- mean(cor2.matrix[targets_buffer,], na.rm = TRUE)
cor2_buffer_sd <- sd(cor2.matrix[targets_buffer,], na.rm = TRUE)
cor2_buffer_cutoff <- cor2_buffer_mean+3*(cor2_buffer_sd)

###Create sample specific buffer means for normalisation
cor2_buffer_sample_mean <- colMeans(cor2.matrix[targets_buffer,], na.rm = TRUE)
log_buffer_sample_mean <- log2(cor2_buffer_sample_mean)

#subtract buffer mean from each sample to generate new intensity matrices (log transformed data)
#remove highly variable samples from analyses
norm.matrix <- log.cor.matrix
for(i in 1:ncol(norm.matrix))
{
  norm.matrix[,i] <- norm.matrix[,i]-log_buffer_sample_mean[i]
}


# Set negative normalized log values to zero...
norm2.matrix <- norm.matrix
i = 1
for (i in 1:length(norm2.matrix))
{
  if (is.na(norm2.matrix[[i]]) | norm2.matrix[[i]] > 0) 
  {
    i = i+1
    } else if(norm2.matrix[[i]] < 0) 
    {
      norm2.matrix[[i]] <- 0
      i = i+1
      }
}

### Average duplicates, if the data has technical replicates in the form of 2 blocks / subarray

if (reps == 2)
{
  n = nrow(norm2.matrix)/2
  rep1 <- norm2.matrix[1:n,]
  rep2 <- norm2.matrix[(n+1):(n*2),]
  
  norm2average.matrix <- matrix(nrow = n, ncol = ncol(norm2.matrix))
  colnames(norm2average.matrix) = colnames(norm2.matrix)
  rownames(norm2average.matrix) = rownames(norm2.matrix[(1:n),])
  
  norm2average.matrix <- (rep1+rep2)/2
  
}
  
  ### Check for deviant technical replicates, automatically exclude (set to NA)
  # Use Patrick’s formula for ELISA to compare replicates within one subarray
  # if rep1 or rep2 is more than 1.5 times rep2 or rep1, respectively, and either is above 0.2, then exclude
  # I have adapted this for the log data, I think it was written for raw MFI.
  # I don't think it makes sense to use 0.2 as a cutoff anymore, I think we should leave in all the data and things will be excluded later as seronegative.
  # Also, can redo this using the subsetted matrices (rep1 and rep2) and it should be shorter
if (reps == 2)
{ 
  for (k in 1:ncol(norm2.matrix))
  {
    for(j in 1:n) 
    {
      if (is.na(norm2.matrix[j,k])|is.na(norm2.matrix[(j+n),k]))
        { 
        j+1
        } else if (norm2.matrix[j,k] > (log2(1.5) + norm2.matrix[(j+n),k]) | (norm2.matrix[(j+n),k] > (log2(1.5) + norm2.matrix[j,k])) == TRUE) 
        {
        norm2average.matrix[j,k] <- NA
        }
    }
  }
}
remove(j,k)

## Calculate correlation coefficient (default is pearson) 
repR <- cor(c(rep1), c(rep2), use = "complete.obs")
print(repR)

## Plot replicate 1 v. replicate 2 for each protein or each person and calculate correlation coefficient, flag r<0.95.
png(filename = paste("replicatescorrelation.tif"), width = 5, height = 4, units = "in", res = 600)
par(mar = c(4, 3, 1, 0.5), oma = c(1, 1, 1, 1), bty = "o", 
    mgp = c(2, 0.5, 0), cex.main = 1, cex.axis = 0.5, cex.lab = 0.7, xpd=NA, las=1)

plot(rep1, rep2, col="red", cex = 0.1)
mtext(c(paste("Pearson correlation coefficient:", round(repR, digits=4))), side=3, adj=0)

graphics.off()

###DATA ANALYSIS###
#Before doing any further analysis, we have to get rid of samples or targets that we are no longer interested in
#E.g. If control individuals are in our analysis, they will affect mixture model based cut-offs
#E.g. If control targets are still in our analysis, they will muck up our protein breadth estimates
#This means we have to subset the data, so some earlier annotations will from here on be wrong (e.g. index_sample will no longer equal 96)
rmsamp_all <- unique(c(targets_blank, targets_buffer, targets_ref, targets_std, high_targets_disinclude))
norm_sub.matrix <- norm.matrix[-rmsamp_all, samples_test]
samples_sub.df <- samples.df[samples_test,]

###Form a seropositivty matrix based on reactivity over the sample background.
#The cut-off is 1. As the data is log2 normalised, a value of 1 equates to eaxctly doible the MFI of the samples buffer mean MFI.
#An alternative would be to raise this threshold to mean + 3SD. I think this is fine.
seropos_buffer_sub.matrix <- as.matrix(norm_sub.matrix > 1)+0

###We may also wish to attempt mixture model cut-offs for each antigen. This is more complex.
#THIS DOESNT WORK
cutoff_list <- c()
restarts_list <-c()

for(i in 1:296){
  temp <-  normalmixEM(t(norm_sub.matrix[i,]), lambda = NULL, mu = NULL, sigma = NULL, maxrestarts=100000)
  if ( temp$mu[1] < temp$mu[2]) cutoff_temp <- temp$mu[1]+3*temp$sigma[1]
  if (temp$mu[1] > temp$mu[2])  cutoff_temp <- temp$mu[2]+3*temp$sigma[2]
  restarts_temp <- temp$restarts
  cutoff_list <- c(cutoff_list,cutoff_temp)
  restarts_list <- c(restarts_list,restarts_temp)
}

seropos_mix.matrix <- norm.matrix
for(i in 1:index_target)
{
  if(restarts_list[i]<=100000){
    seropos_mix[,i] <- as.numeric((norm.matrix.df[i,] > cutoff_list[i]))
  } else {
    seropos_mix[,i]<-rep(NA,dim(seropos_mix)[1])
  }  
}

###Create a threshold for overall target and person reactivity
#e.g. To be included in heatmaps and other analyses, perhaps targets should be reacted to by at least 5% of people?
#Similarly, perhaps unreactive individuals should be disincluded? Either way - this is informative.
target_breadth <- rowSums(seropos_buffer_sub.matrix, na.rm=TRUE)
target_reactive <- target_breadth > (ncol(seropos_buffer_sub.matrix)/100)*5
cat(sum(target_reactive), "out of", nrow(seropos_buffer_sub.matrix), "protein targets are reactive in at least 5% of people")

person_breadth <- colSums(seropos_buffer_sub.matrix, na.rm=TRUE)
person_exposed <- person_breadth > (nrow(seropos_buffer_sub.matrix)/100)*5
cat(sum(person_exposed), "out of", ncol(seropos_buffer_sub.matrix), "samples are reactive to at least 5% of proteins")

###Create person and protein magnitudes
#Not subset - i.e. taking the mean magnitude of response to all 296 target proteins
png(filename = paste("magnitude_breadth_chmitimepoints.tif"), width = 6, height = 3, units = "in", res = 600)
par(mfrow=c(1,3), mar = c(4, 3, 1, 0.5), oma = c(1, 1, 1, 1), bty = "o", 
    mgp = c(2, 0.5, 0), cex.main = 1, cex.axis = 0.5, cex.lab = 0.7, xpd=NA, las=1)
person_magnitude1 <- colMeans(norm_sub.matrix[,])
boxplot(person_magnitude1~samples_sub.df$Timepoint, ylab="Mean magnitude of response to all proteins", xlab="Timepoint", col="transparent")
#Subset only using reactive targets (60/296)
person_magnitude2 <- colMeans(norm_sub.matrix[target_reactive==TRUE,])
boxplot(person_magnitude2~samples_sub.df$Timepoint, ylab="Mean magnitude of response to all proteins", xlab="Timepoint", col="transparent")
###Some boxplots based on the available data
samples_sub.df$Timepoint<-factor(samples_sub.df$Timepoint, levels=c("d0", "C-1", "C+28", "C+35"))
boxplot(person_breadth~samples_sub.df$Timepoint, ylab="Breadth of response", xlab="Timepoint", col="transparent")
graphics.off()

###Subset the final corrected dataframe to only look at relevant target spots
rmsamp_all <- unique(c(targets_blank, targets_buffer, targets_ref, targets_std, high_targets_disinclude))
cor4.matrix <- cor3.matrix[-rmsamp_all,]

###Attemp at plot showing treatment stages
png(filename = paste("baseline_correct4.tif"), width = 8, height = 3, units = "in", res = 600)
par(mfrow=c(1,5), mar = c(4, 3, 1, 0.5), oma = c(1, 1, 1, 1), bty = "o", 
    mgp = c(2, 0.5, 0), cex.main = 1, cex.axis = 0.5, cex.lab = 0.7, xpd=NA, las=1)
boxplot(rowMeans(fore.matrix[-targets_allcontrol,]), fg="blue")
boxplot(rowMeans(cor.matrix[-targets_allcontrol,]))
boxplot(rowMeans(cor2.matrix[-targets_allcontrol,]))
boxplot(rowMeans(cor3.matrix[-targets_allcontrol,]))
boxplot(rowMeans(norm.matrix[-targets_allcontrol,]))
graphics.off()

#Plot the corrected values - comparing between correction methods
png(filename = paste("baseline_correct.tif"), width = 5, height = 4, units = "in", res = 600)
par(mfrow=c(3,1), mar = c(4, 3, 1, 0.5), oma = c(1, 1, 1, 1), bty = "o", 
    mgp = c(2, 0.5, 0), cex.main = 1, cex.axis = 0.5, cex.lab = 0.7, xpd=NA, las=1)

hist(cor.matrix, add=FALSE, main=paste("Corrected MFI"), breaks=100, xlab=NA)
hist(cor2.matrix, add=FALSE, main=paste("Blunt corrected MFI"), breaks=100, xlab=NA)
hist(cor3.matrix, add=FALSE, main=paste("LIMMA corrected MFI"), breaks=100, xlab="MFI")
graphics.off()

png(filename = paste("baseline_correct2.tif"), width = 5, height = 4, units = "in", res = 600)
par(mar = c(4, 3, 1, 0.5), oma = c(1, 1, 1, 1), bty = "o", 
    mgp = c(2, 0.5, 0), cex.main = 1, cex.axis = 0.5, cex.lab = 0.7, xpd=NA, las=1)
plot(log(cor2.matrix), log(cor3.matrix))
graphics.off()

png(filename = paste("baseline_correct5.tif"), width = 5, height = 4, units = "in", res = 600)
par(mfrow=c(1,5), mar = c(4, 3, 1, 0.5), oma = c(1, 1, 1, 1), bty = "o", 
    mgp = c(2, 0.5, 0), cex.main = 1, cex.axis = 0.5, cex.lab = 0.7, xpd=NA, las=1)
a <- density(fore.matrix)
b <- density(cor.matrix)
c <- density(cor2.matrix)
d <- density(cor3.matrix)
e <- density(norm.matrix, na,rm=TRUE)
plot.density(a, main=paste("Foreground MFI"))
plot.density(b, main=paste("Background corrected MFI"))
plot.density(c, main=paste("Baseline corrected MFI 1"))
plot.density(d, main=paste("Baseline corrected MFI 2"))
plot.density(e, main=paste("Normalised MFI"))
graphics.off()

graphics.off()


#Test heatmap
#THIS DOESN'T WORK

# install.packages("pheatmap", "RColorBrewer", "viridis")
library(pheatmap)
library(RColorBrewer)
library(viridis)

# Data frame with column annotations.
mat_col <- data.frame(group=samples_sub.df$Timepoint)
rownames(mat_col) <- colnames(norm_sub.matrix)

# List with colors for each annotation.
mat_colors <- list(group = brewer.pal(4, "Set1"))
names(mat_colors$group) <- unique(samples_sub.df$Timepoint)

pheatmap(
  mat               = norm_sub.matrix,
  color             = inferno(10),
  border_color      = NA,
  show_colnames     = FALSE,
  show_rownames     = FALSE,
  annotation_col    = mat_col,
  annotation_colors = mat_colors,
  drop_levels       = TRUE,
  fontsize          = 14,
  main              = "Default Heatmap"
)