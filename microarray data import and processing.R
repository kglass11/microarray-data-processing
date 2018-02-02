###Combined script for reading in and processing microarray data to prepare for analysis

###Create a folder, into which you should copy this script, your .gpr files (named 'Slide 1.gpr', 'Slide 2.gpr' etc.), and your sample list csv file.
# Your sample list file needs four columns. Their names must be exactly as written here, though the order of the samples does not matter:
#1.slide_no
#2.sample_id
#3.block_rep_1
#4.block_rep_2

###Clear the environnment - OR go to Session > Clear Workspace
rm(list=ls())

###Install any packages you may need for this script. Go to Tools, install packages, 
#then install from CRAN repository or local file or run:
#install.packages(c("gtools","contrast", "beeswarm", "mixtools", "gplots", "ggplot2", "gcookbook"))

##Install limma if you haven't already
## try http:// if https:// URLs are not supported
#source("https://bioconductor.org/biocLite.R")
#biocLite("limma")

###Load packages needed for this script
require("gtools")

library(limma)
library(contrast)
library(beeswarm)
library(mixtools)
library(gplots)
library(ggplot2)
library(gcookbook)

### Define variables based on your study that will be used later in the script
# define working directory character vector, example "I:/Drakeley Group/Protein microarrays/Experiments/310817 Bijagos Islands/Screen 2"
workdir <- "/Users/Katie/Desktop/R files from work/PRISM Immune v Nonimmune"

# define a shorthand name for your study which will be appended in the file name of all exported files
study <- "PRISM1"

#define file name for sample IDs character vector, example "Analysis sample list 2.csv"
sample_file <- "Analysis sample list.csv"

#define file name for sample file + additional metadata 
#(ex: test v. control sample, exclude samples (based on slide image), experimental groups)
meta_file <- "Immune v nonimmune metadata.csv"

#number of technical replicates for the study (usually 1 or 2)
reps <- 2

#define number of blocks per slide
index_block <- 32

### Set the working directory for the folder where .gpr files are. Can check working
#directory after with getwd()
setwd(workdir)

#####################################
###READING IN YOUR MICROARRAY data###
#####################################

###Identify all .gpr files in your set path. default path for list.files() is the working directory.
slide_ids <- list.files(pattern="*.gpr")

###Read in the contents of each of the .gpr files you just identified
#This simply combines the data in a list format
#If your gpr files differ in terms of the number of rows of information before the header, you need to edit this command
#i.e skip=34 means it will skip the first 34 rows, and header=T means it will read the 35th ow as a header, and the 36th as the first row of data
# which is what currently happens in our GPR files.
#Note: slide_no_temp should be the last slide added to the list
#Slide number is assigned based on the file name. As files are created without editing by the genepix software - this seems like the easiest way of identifying files.
slides_list <- list()
for(i in 1:length(slide_ids)) { 
  
  slides_list[[i]] <- read.table(slide_ids[i],skip=34,sep="\t",header=T)
  slide_no_temp <- substr(slide_ids[[i]],7,nchar(slide_ids[[i]])-4)
  slides_list[[i]][,"slide_no"] <- slide_no_temp
}
remove(i)

###Name your slides in this list, according to their file names
names(slides_list) <- slide_ids

###Bind all data from the slide data list (slides.list) into a single dataframe
slides_all.df <- c()
for(i in 1:length(slides_list)) { 
  
  slides_all.df <- rbind(slides_all.df,slides_list[[i]])
}
remove(i)

###Read in list of sample IDs
samples.df <- read.csv(sample_file, header=T, na.strings = " ", check.names = FALSE, stringsAsFactors = FALSE)

###Create a vector listing all of your samples, in the order they appear in your samples_list file
#In our samples_list files, the sample ID is always in column two - which is why that column is picked out here
samples <- as.character(samples.df[1:nrow(samples.df), 2])

###Now, make a new sample variable that is unique (to avoid issues with multiple blanks, etc.) 
# in the next version of Will's, he changed the samples_unique to column 11 instead of 5... why? 
# we need this not to change based on the sample file/study, it should be the same every time
samples_unique <- c(paste(as.character(samples.df[1:nrow(samples.df), 2]), rownames(samples.df), sep = "_"))
samples.df[,5] <- samples_unique
colnames(samples.df)[5] <- "sample_id_unique"

###Create vectors indicating the number of slides, blocks, and samples
#Slide and sample number are determined automatically from the data you input, whereas block number is manual in this instance
index_slide <- as.numeric(length(slides_list))
index_sample <- as.numeric(length(samples))

###Assign your sample_ids to each row of the combined slide data (slides_all.df)
#The order of data in your samples.df file is irrelevant, as long as each sample ID is correctly matched to its slide and block numbers
slides_all.df$sample_id_unique <- c()

for(i in 1:dim(slides_all.df)[1]){
  
  print(i)
  row_ite<-slides_all.df[i,]
  block_ite<-row_ite$Block
  slide_ite<-row_ite$slide_no
  sample_info_1<-samples.df[which(samples.df$slide_no==slide_ite),]
  match<-block_ite%in%sample_info_1[,"block_rep_1"]
  if(match==TRUE){
    value<-which(block_ite==sample_info_1[,"block_rep_1"])
  }else{
    value<-which(block_ite==sample_info_1[,"block_rep_2"])
  }
  sample_info_2<-sample_info_1[value,"sample_id_unique"]
  slides_all.df$Sample[i]<-as.character(sample_info_2)
}
remove(i, sample_info_1, sample_info_2, match, row_ite, block_ite, slide_ite)

###Write slides_all.df to a file to keep as a csv in your directory
write.csv(slides_all.df,file=paste0(study,"_slidesall_combinedGPR.csv"), row.names=T)

### Make a spot annotations dataframe
#Column 42 is the sample id column
annotation_targets.df <- slides_all.df[which(slides_all.df[,42]==samples_unique[1]),1:4]
annotation_targets.df <- cbind(row.names(annotation_targets.df), annotation_targets.df)
colnames(annotation_targets.df)[1] <- "target_id_numeric"
annotation_targets.df[6] <- c(paste(rownames(annotation_targets.df), annotation_targets.df$Name, annotation_targets.df$Block, sep = "_"))
colnames(annotation_targets.df)[6] <- "target_id_unique"
rownames(annotation_targets.df) <- c(annotation_targets.df$target_id_unique)

###Make a final index, this time of targets
index_target <- as.numeric(length(annotation_targets.df$target_id_numeric))

#####################################
###PROCESSING YOUR MICROARRAY data###
#####################################

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
write.csv(cor.matrix, file=paste0(study,"_background_corrected.matrix.csv"))

###Assign target names to groups of your array targets to identify their 'type'
targets_blank = c(grep("BLANK", annotation_targets.df$Name))
targets_buffer = c(grep("buffer", annotation_targets.df$Name))
targets_ref = c(grep("REF", annotation_targets.df$Name))
targets_std = c(grep("Std", annotation_targets.df$Name))
targets_allcontrol = c(targets_blank, targets_buffer, targets_ref, targets_std)

###Remove "bad" spots from subsequent analysis

#First, we need to tell R which spots are the highest spots on the plate
#the spots printed immediately after these can be made to have null MFIs
#This call will create a vector identifying all ref spots, plus all std spots called "Std 1"
high_targets <- c(targets_ref, grep("Std 1", row.names(annotation_targets.df)))

#To identify the spots to disinclude, we need to identify the position of the next spot in the print run.
#Because block 2 is printed before block 1, we can subtract an entire block, minus 1 row (e.g. 180-12=168)
high_targets_disinclude <- ifelse(high_targets>=index_target/2, high_targets-((index_target/2)-12), high_targets+index_target)
high_targets_disinclude <- high_targets_disinclude[which(high_targets_disinclude<=(index_target-24))] 

###Convert the spots to be disincluded to NAs in background corrected data
cor2.matrix <- cor.matrix
cor2.matrix[high_targets_disinclude, ] <- NA

###QUALITY CONTROL###

### 1. Background
#To identify which slides, pads, and samples are significantly deviated, we need to calculate the mean and generate an arbitrary cut off 
#The cut off can be used to flag samples, and tells us whether deviation is universal, or specific to slides, pads, or samples
#We are not automatically excluding samples with deviant background, but it might be a good idea to
# go back and review the images for those pads.

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
write.csv(back_sample_deviant, file = paste0(study,"_deviant_sample_background.csv"))
#Coefficient of variation for all background datapoints, and all exlcuding deviant samples
back_cov_all <- sd(back.matrix)/mean(back.matrix)
back_cov_normal <- sd(back.matrix[, back_normal])/mean(back.matrix[, back_normal])

#Plots by slide/pad/sample
png(filename = paste0(study,"_background_mfi_samples.tif"), width = 5.5, height = 10, units = "in", res = 600)
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
png(filename = paste0(study,"_background_mfi_targets.tif"), width = 15, height = 10, units = "in", res = 600)
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
png(filename = paste0(study,"_cov_back_sample.tif"), width = 10, height = 4, units = "in", res = 600)
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
#The mean and cutoff will be calculated EXCLUDING any "bad" buffer targets previously set to NA.
#The cut off can be used to flag samples, and tells us whether deviation is universal, or specific to slides, pads, or samples

#***Samples which have deviant buffer means will be automatically excluded from further analysis
#***Buffer targets which are deviant across all pads and slides will be automatically excluded from further analysis

#Buffer assessment EXCLUDING "bad" spots: 
#Mean median corrected MFI for all data points
cor_buffer_mean <- mean(cor2.matrix[targets_buffer,], na.rm = TRUE)
#SD median corrected MFI for all data points
cor_buffer_sd <- sd(cor2.matrix[targets_buffer,], na.rm = TRUE)
#Cut-off for deviation from mean
cor_cutoff <- cor_buffer_mean+(3*cor_buffer_sd)

###Identify deviant samples/slides/pads

# EXCLUDING "bad" spots - Generate a vector of mean buffer intensity for EACH sample (corrected person magnitude) 
cor_buffer_sample_mean <- colMeans(cor2.matrix[targets_buffer,], na.rm = TRUE)
#Vectors to identify normal and deviant samples
cor_normal <- which(cor_buffer_sample_mean<=cor_cutoff)
cor_deviant <- which(cor_buffer_sample_mean>cor_cutoff)
#Generate a table showing which slides, pads, samples have deviant corrected buffer values
cor_sample_deviant <- samples.df[cor_deviant,]
cor_sample_deviant
write.csv(cor_sample_deviant, file = paste0(study,"_deviant_sample_buffer.csv"))

#Coefficient of variation for buffer datapoints (EXCLUDING "bad" spots), and exlcuding deviant samples
cor_buffer_cov <- cor_buffer_sd/cor_buffer_mean
cor_buffer_cov_normal <- sd(cor2.matrix[targets_buffer, cor_normal], na.rm = TRUE)/mean(cor2.matrix[targets_buffer, cor_normal], na.rm = TRUE)

#Plot CoV of buffer by sample - EXCLUDING "bad" spots
cor_buffer_cov_sample <- c()

png(filename = paste0(study,"_cov_buffer.tif"), width = 5, height = 4, units = "in", res = 600)
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

#Buffer assessment INCLUDING "bad" spots:
cor_buffer_all_mean <- mean(cor.matrix[targets_buffer,])
#SD median corrected MFI for all data points
cor_buffer_all_sd <- sd(cor.matrix[targets_buffer,])
#Cut-off for deviation from mean
cor_all_cutoff <- cor_buffer_all_mean+(3*cor_buffer_all_sd)

#Coefficient of variation for all buffer datapoints (INCLUDING "bad" spots), and all exlcuding deviant samples
cor_buffer_cov_all <- cor_buffer_all_sd/cor_buffer_all_mean
cor_buffer_cov_all_normal <- sd(cor.matrix[targets_buffer, cor_normal])/mean(cor.matrix[targets_buffer, cor_normal])

#Plots by slide/pad/sample INCLUDING "bad" spots (all buffer data)
png(filename = paste0(study,"_buffer_mfi_QCplots.tif"), width = 5.5, height = 10, units = "in", res = 600)
par(mfrow=c(3,1), mar = c(2, 3, 2.25, 0.5), oma = c(11.5, 0, 1, 0), bty = "o", 
    mgp = c(2, 0.5, 0), cex.main = 1.5, cex.axis = 0.6, cex.lab = 0.9, xpd=NA, las=1)
boxplot(t(cor.matrix[targets_buffer,]) ~ samples.df$slide_no, outcex=0.5,
        ylab="Corrected MFI", xlab="Slide", add=FALSE)
abline(h = cor_all_cutoff, col = "red", lty = 2, lwd = 0.7, xpd=FALSE)
title(main = "Corrected MFI (buffer only) by SLIDE/SUBARRAY/SAMPLE\n", adj=0)
boxplot(t(cor.matrix[targets_buffer,]) ~ samples.df$block_rep_1, outcex=0.5,
        ylab="Corrected MFI", xlab="Subarrays (1/2, 3/4, etc.)", add=FALSE, las=1)
abline(h = cor_all_cutoff, col = "red", lty = 2, lwd = 0.7, xpd=FALSE)
boxplot(t(cor.matrix[targets_buffer,]) ~ samples.df$sample_id_unique, outcex=0.5,
        ylab="Corrected MFI", xlab="Sample", add=FALSE, las=2, cex.axis = 0.4, yaxt="n")
axis(2, cex.axis=0.6)
abline(h = cor_all_cutoff, col = "red", lty = 2, lwd = 0.7, xpd=FALSE)

mtext(c(paste("INCLUDING Bad Buffer Spots / EXCLUDING Bad Buffer Spots:")), side=1, cex=0.8, line=2, outer=TRUE, xpd=NA, adj=0)
mtext(c(paste("Mean overall corrected buffer MFI:", round(cor_buffer_all_mean, digits=3), "/", round(cor_buffer_mean, digits=3))), side=1, cex=0.8, line=3.5, outer=TRUE, xpd=NA, adj=0)
mtext(c(paste("SD overall corrected buffer MFI:", round(cor_buffer_all_sd, digits=3), "/", round(cor_buffer_sd, digits=3))), side=1, cex=0.8, line=5, outer=TRUE, xpd=NA, adj=0)
mtext(c(paste("Cut-off overall corrected buffer MFI:", round(cor_all_cutoff, digits=3), "/", round(cor_cutoff, digits=3))), side=1, cex=0.8, line=6.5, outer=TRUE, xpd=NA, adj=0)
mtext(c(paste("CoV overall corrected buffer MFI:", round(cor_buffer_cov_all, digits=3), "/", round(cor_buffer_cov, digits=3))), side=1, cex=0.8, line=8, outer=TRUE, xpd=NA, adj=0)
mtext(c(paste("CoV overall corrected buffer MFI (excl. deviant samples):", round(cor_buffer_cov_all_normal, digits=3), "/", round(cor_buffer_cov_normal, digits=3))), side=1, cex=0.8, line=9.5, outer=TRUE, xpd=NA, adj=0)

graphics.off()

###Identify deviant buffer targets across all pads

#EXCLUDING "bad" spots - Mean corrected MFI for every sample for EACH target (background protein magnitude)
cor_target_mean <- rowMeans(cor2.matrix)
#Are any corrected buffer targets universally deviant, accross all pads?
deviant_buffer_targets <- Reduce(intersect, list(targets_buffer, which(cor_target_mean>cor_cutoff)))
cor_buffer_deviant.df <- annotation_targets.df[deviant_buffer_targets,]
cor_buffer_deviant.df
write.csv(cor_buffer_deviant.df, file = paste0(study,"_deviant_buffer_targets.csv"))

#list of disincluded buffer targets ("bad" spots)
removed_buffer_targets <- c()
for (i in 1:length(targets_buffer)){
  for(j in 1:length(high_targets_disinclude))
    if (targets_buffer[i] == high_targets_disinclude[j]){
      removed_buffer_targets <- targets_buffer[i]
    }
}
remove(i,j)

#All slides, INCLUDING "bad" spots (ALL buffer targets)
png(filename = paste0(study, "_buffer_targets.tif"), width = 5, height = 4, units = "in", res = 600)
par(mar = c(5, 3, 2.25, 0.5), oma = c(0, 0, 0, 0), bty = "o", 
    mgp = c(2, 0.5, 0), cex.main = 1, cex.axis = 0.6, cex.lab = 0.7, xpd=NA, las=2)
boxplot(t(cor.matrix[targets_buffer,]),
        cex=0.5,
        ylab="Corrected MFI")
abline(h = cor_cutoff, col = "red", lty = 2, lwd = 0.7, xpd=FALSE)

mtext(c(paste("Excluded buffer targets:", removed_buffer_targets)), las = 1)

graphics.off()

#By slide for up to 6 slides (the last slides) - INCLUDING "bad" spots (ALL buffer targets)
#mfrow is ordered by the number of rows, and columns of plots you want - so must be edited based on the number of slides
png(filename = paste0(study, "_buffer_spots_slide.tif"), width = 5, height = 4, units = "in", res = 600)
par(mfrow=c(2,3), mar = c(4, 3, 1, 0.5), oma = c(1, 1, 1, 1), bty = "o", 
    mgp = c(2, 0.5, 0), cex.main = 1, cex.axis = 0.5, cex.lab = 0.7, xpd=NA, las=2)
for (i in 1:index_slide){
  boxplot(t(cor.matrix[targets_buffer,samples.df$slide_no==i]),
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

#Automatically exclude the samples with deviant buffer values in a new matrix
cor3.matrix <- cor2.matrix
cor3.matrix[,cor_deviant] <- NA

#Automatically exclude the buffer targets deviant across all arrays
#KG - It would be good to test this on data that has deviant buffer targets! 
cor3.matrix[deviant_buffer_targets,] <- NA

### Mean values for each target ordered by position within the arrays (not log-transformed or normalized data)
#KG - I think we might want this somehwere else in the script using a more processed matrix

#Check average corrected values for each target for all individuals
cor_target_mean <- rowMeans(cor3.matrix, na.rm = TRUE)
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
write.csv(cbind(cor_target_mean_b1b2, annotation_target_b1b2), file=paste0(study,"_target_mean_as_array.csv"))
remove(cor_target_mean, cor_target_mean_b1, cor_target_mean_b2,cor_target_mean_b1b2, annotation_target_b1, annotation_target_b2, annotation_target_b1b2)


#### LOG TRANSFORMATION AND NORMALIZATION ###

### Log transform the data (base 2)
log.cor.matrix <- log2(cor3.matrix)

### Normalization

###Create sample specific buffer means for normalisation
cor2_buffer_sample_mean <- colMeans(cor3.matrix[targets_buffer,], na.rm = TRUE)
cor2_buffer_sample_sd <- apply(cor3.matrix[targets_buffer,], 2, sd, na.rm = TRUE)
#log2 transform sample buffer mean
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

write.csv(norm2.matrix, file = paste0(study,"_normalized_log_data.csv"))
