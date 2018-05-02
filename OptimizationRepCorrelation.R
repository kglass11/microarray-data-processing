#continuing on from optimizationanalysis.R and OptimizationTop10.R

#Replicates correlation and other correlation analysis

#1. Correlation replicate 1 with replicate 2 for all conditions (480 conditions) 

  #rep1.matrix and rep2.matrix have the data with negative normalized values set to zero
  #trans.norm.rep1 and trans.norm.rep2 have the data without setting negative normalized values to zero
  #trans.norm.meta.rep1 and 2 have the data with sample_id as a variable not rownames and have 
    #slide_type, blocking_buffer, and block_dilution as final 3 columns
  #data is not a ratio of positive to negative, nor has GST been subtracted
  #data includes buffer and ref spots as well

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

Rep2.df$sample <- "CP3"
Rep2.df$sample[c(grep("PRISM", Rep2.df$sample_id))] <- "PRISM"
Rep2.df$sample[c(grep("Swazi", Rep2.df$sample_id))] <- "Swazi"

#the names of the columns have spaces. So they cannot be used in lmer and other functions
newnames <- make.names(colnames(Rep1.df))
colnames(Rep1.df) <- newnames

newnames <- make.names(colnames(Rep2.df))
colnames(Rep2.df) <- newnames


  #calculate pearson's correlation coefficients for replicate 1 vs replicate 2 row-wise

