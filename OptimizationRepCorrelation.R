#continuing on from optimizationanalysis.R and OptimizationTop10.R

#saved after corrgrams
save.image("RepsCor.RData")

load("RepsCor.RData")

#Replicates correlation and other correlation analysis

# install.packages("corrgram")
library(corrgram)

###1. Correlation replicate 1 with replicate 2 for all conditions (480 conditions) 

  #rep1.matrix and rep2.matrix have the data with negative normalized values set to zero
  #trans.norm.rep1 and trans.norm.rep2 have the data without setting negative normalized values to zero
  #trans.norm.meta.rep1 and 2 have the data with sample_id as a variable not rownames and have 
    #slide_type, blocking_buffer, and block_dilution as final 3 columns
  #data is not a ratio of positive to negative, nor has GST been subtracted
  #data includes buffer and ref spots as well, but exclude these when separating by print buffer.
    #We don't want to count those because they are printed in system buffer only, and don't have a specific print buffer

  #Separate data by print buffers for replicate 1 and replicate 2.

#edit the below script for trans.norm.rep1 and trans.norm.rep2

#make a target_meta.df with the block 2 antigen names
target_meta2.df <- target_meta.df
target_meta2.df$Name <- colnames(trans.norm.rep2)

#merge transposed final data with target metadata again
target1.df <- merge(target_meta.df, t(trans.norm.rep1), by.y="row.names", by.x = "Name", sort = FALSE)
target2.df <- merge(target_meta2.df, t(trans.norm.rep2), by.y="row.names", by.x = "Name", sort = FALSE)

#separate data by print buffer in three data frames
PB1.1 <- filter(target1.df, Print_Buffer == "1")
PB2.1 <- filter(target1.df, Print_Buffer == "2")
PB3.1 <- filter(target1.df, Print_Buffer == "3")

PB1.2 <- filter(target2.df, Print_Buffer == "1")
PB2.2 <- filter(target2.df, Print_Buffer == "2")
PB3.2 <- filter(target2.df, Print_Buffer == "3")

#replace row names with new rownames that are the same for all print buffers, don't include spaces
rownames(PB1.1) <- PB1.1$Name
rownames(PB2.1) <- PB1.1$Name
rownames(PB3.1) <- PB1.1$Name

rownames(PB1.2) <- PB1.2$Name
rownames(PB2.2) <- PB1.2$Name
rownames(PB3.2) <- PB1.2$Name

#Print buffers: 1 = AJ Buffer C		2 = AJ Glycerol buffer		3 = Nexterion Spot
PB1.1 <- as.data.frame(t(PB1.1[,6:ncol(PB1.1)]))
PB1.1$print_buffer <- "AJ"
PB2.1 <- as.data.frame(t(PB2.1[,6:ncol(PB2.1)]))
PB2.1$print_buffer <- "AJ_Glycerol"
PB3.1 <- as.data.frame(t(PB3.1[,6:ncol(PB3.1)]))
PB3.1$print_buffer <- "Nexterion_Spot"

PB1.2 <- as.data.frame(t(PB1.2[,6:ncol(PB1.2)]))
PB1.2$print_buffer <- "AJ"
PB2.2 <- as.data.frame(t(PB2.2[,6:ncol(PB2.2)]))
PB2.2$print_buffer <- "AJ_Glycerol"
PB3.2 <- as.data.frame(t(PB3.2[,6:ncol(PB3.2)]))
PB3.2$print_buffer <- "Nexterion_Spot"

#merge these with sample_meta_f to get other columns, then rbind all together
PB1.1.meta <- merge(sample_meta_f.df, PB1.1, by.x = "sample_id", by.y = "row.names",sort=FALSE)
PB2.1.meta <- merge(sample_meta_f.df, PB2.1, by.x = "sample_id", by.y = "row.names",sort=FALSE)
PB3.1.meta <- merge(sample_meta_f.df, PB3.1, by.x = "sample_id", by.y = "row.names",sort=FALSE)

PB1.2.meta <- merge(sample_meta_f.df, PB1.2, by.x = "sample_id", by.y = "row.names",sort=FALSE)
PB2.2.meta <- merge(sample_meta_f.df, PB2.2, by.x = "sample_id", by.y = "row.names",sort=FALSE)
PB3.2.meta <- merge(sample_meta_f.df, PB3.2, by.x = "sample_id", by.y = "row.names",sort=FALSE)

Rep1.df <- rbind(PB1.1.meta, PB2.1.meta, PB3.1.meta)
Rep2.df <- rbind(PB1.2.meta, PB2.2.meta, PB3.2.meta)

#need to change the sample_id to only CP3, PRISM, and Swazi 
#so that the model knows they are the same sample
Rep1.df$sample <- "CP3"
Rep1.df$sample[c(grep("PRISM", Rep1.df$sample_id))] <- "PRISM"
Rep1.df$sample[c(grep("Swazi", Rep1.df$sample_id))] <- "Swazi"
Rep1.df$sample[c(grep("Neg", Rep1.df$sample_id))] <- "Neg"

Rep2.df$sample <- "CP3"
Rep2.df$sample[c(grep("PRISM", Rep2.df$sample_id))] <- "PRISM"
Rep2.df$sample[c(grep("Swazi", Rep2.df$sample_id))] <- "Swazi"
Rep2.df$sample[c(grep("Neg", Rep2.df$sample_id))] <- "Neg"

#the names of the columns have spaces. So they cannot be used in lmer and other functions
newnames <- make.names(colnames(Rep1.df))
colnames(Rep1.df) <- newnames

newnames <- make.names(colnames(Rep2.df))
colnames(Rep2.df) <- newnames

#calculate pearson's correlation coefficients for replicate 1 vs replicate 2 row-wise
#there are 30 rows for which all columns have NA, so correlation matrix also has NA for those conditions
repcor <- c(nrow(Rep1.df))

for(i in 1:nrow(Rep1.df)){
  repcor[i] <- cor(x = c(as.matrix(Rep1.df[i,11:46])), y = c(as.matrix(Rep2.df[i,11:46])), use = "everything")
}
remove(i)

#add correlation coefficients to the metadata, but have correlation coeffs for Neg data as well
repcor.df <- as.data.frame(repcor)
repcor.df$sample_id <- Rep1.df$sample_id
repcor.meta.df <- merge(sample_meta_f.df, repcor.df, sort = TRUE)

#add print buffer and sample - *** need to make sure these are the right print buffer and sample!!
sortRep1 <- Rep1.df[order(Rep1.df$sample_id),]
repcor.meta.df$print_buffer <- sortRep1$print_buffer
repcor.meta.df$sample <- sortRep1$sample

#set factor order for sample
repcor.meta.df$sample <- factor(repcor.meta.df$sample, levels = as.character(c("CP3", "PRISM", "Swazi", "Neg")))

#plot All correlation coefficients vs Slide number
png(filename = paste0("RepCor1.tif"), width = 8, height = 3.5, units = "in", res = 1200)
par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))

ggplot(repcor.meta.df, aes(x = sample_id, y = repcor))  + 
  geom_point(color = "blue", size = 0.5) + theme_bw() + theme(axis.text.x=element_blank(), axis.ticks.x = element_blank()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(x = "Slide", y = "Correlation Coefficient", title = "Correlation between replicates for each condition and sample")

graphics.off()  

#plot correlation coefficient vs. slide type
png(filename = paste0("RepCor.ST.tif"), width = 8, height = 3.5, units = "in", res = 1200)
par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))

ggplot(repcor.meta.df, aes(x = slide_type, y = repcor, color = slide_type))  + 
  geom_jitter(size = 0.5) + theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 9)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(x = "Slide Type", y = "Correlation Coefficient", title = "Correlation between replicates")

graphics.off()  

#plot correlation coefficient vs. slide type and print_buffer
png(filename = paste0("RepCor.ST.PB.tif"), width = 8, height = 3.5, units = "in", res = 1200)
par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))

ggplot(repcor.meta.df, aes(x = slide_type, y = repcor, color = print_buffer))  + 
  geom_jitter(size = 0.5) + theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 9)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(x = "Slide Type", y = "Correlation Coefficient", title = "Correlation between replicates") +
  geom_violin(aes(x = slide_type, y=repcor))

graphics.off()  

#plot correlation coefficient vs. slide type and sample
png(filename = paste0("RepCor.ST.Sample.tif"), width = 8, height = 3.5, units = "in", res = 1200)
par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))

ggplot(repcor.meta.df, aes(x = slide_type, y = repcor, color = sample))  + 
  geom_jitter(size = 0.5) + theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 9)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(x = "Slide Type", y = "Correlation Coefficient", title = "Correlation between replicates") +
  geom_violin(aes(x = slide_type, y=repcor))

graphics.off() 

#Plot correlation coefficient between replicates for "best" conditions for CP3 only
  #CP3alone4 has the 8 common Top 10% conditions

rep8 <- merge(CP3alone4, repcor.meta.df)
rep8melt <- melt(rep8)
rep8sub <- filter(rep8melt, variable == "repcor")

rep8sub$rowid <- factor(rep8sub$rowid, levels = rev(unique(Top8rowidorder)))

png(filename = paste0("RepCor.Common8.tif"), width = 4.2, height = 2.7, units = "in", res = 1200)

ggplot(rep8sub, aes(x = rowid, y = value, color = sample)) + theme_bw() + 
  geom_point(shape = 18, size = 2) + theme(axis.text.y = element_text(size = 10)) +
  labs(x = "Condition", y = "Correlation Coefficient") +
  scale_color_hue(name = "Sample") +
  theme(axis.text.x = element_text(color = "black"), panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank()) + 
  coord_flip()

graphics.off()

#calculate how many arrays have correlation coefficients above certain thresholds

#function to calculate the percent of correlation coefficients above a certain threshold
repcor.per <- function(threshold, vector){
  pepper <- (vector > threshold) + 0
  x <- sum(pepper, na.rm=TRUE)/length(vector)*100
  return(paste(round(x,1), "% of correlations coefficients are above", threshold))
}

repcor.per(0.95, repcor)
repcor.per(0.9, repcor)
repcor.per(0.8, repcor)
repcor.per(0.7, repcor)
repcor.per(0.5, repcor)

#repeat calculations but for CP3 only
CP3repcor.df <- filter(repcor.meta.df, sample == "CP3")
CP3repcor <- c(CP3repcor.df$repcor)

repcor.per(0.95, CP3repcor)
repcor.per(0.9, CP3repcor)
repcor.per(0.8, CP3repcor)
repcor.per(0.7, CP3repcor)
repcor.per(0.5, CP3repcor)

#Overall correlations between different conditions for the same points
#CP3 "Best" conditions for 4 antigens only - correlogram
#data frame hello has all the right data, but may need to transpose and subset

hellocor <- as.data.frame(t(hello[,13:48]))
colnames(hellocor) <- hello$rowid

png(filename = paste0("Common8.Correlogram.tif"), width = 3.5, height = 3.5, units = "in", res = 1200)
par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))

corrgram(hellocor, order=NULL, lower.panel=panel.pie,
         upper.panel=panel.pts, text.panel=panel.txt)

graphics.off()
