#Sanger Data Analysis Script
#Katie Glass

#####################################
############DATA ANALYSIS############
#####################################

#If you haven't just continued from the processing script, run:
rm(list=ls())

#I:/Drakeley Group/Protein microarrays/Experiments/100817 Sanger/Sanger Data Processed
setwd("/Users/Katie/Desktop/R files from work/100817 Sanger/Sanger Data Processing KG")
getwd()

require("gtools")

library(limma)
library(contrast)
library(beeswarm)
library(mixtools)
library(gplots)
library(ggplot2)
library(gcookbook)
library(dplyr)

#load(file="SangerAfterProcessing.RData")

#install.packages("reshape2")

library(reshape2)

load(file="Sanger.2.Update.RData")

###Seropositivity and Reactivity Thresholds###

#For this section, using normalized data with negative values set to 0 (norm4.matrix)

#Before doing any further analysis, we have to get rid of samples or targets that we are no longer 
#interested in.
#E.g. If control individuals are in our analysis, they will affect mixture model based cut-offs
#E.g. If control targets are still in our analysis, they will muck up our protein breadth estimates
#KG - for now this still includes the same protein target at different dilutions
#This means we have to subset the data, so some earlier annotations will from here on be wrong 
#(e.g. index_sample will no longer equal 96)

#At this point, Remove control samples for further analysis
norm_sub6.df <- norm_sub4.df[,colnames(norm_sub4.df) %in% samples_test]

#Make the dilution column of target_meta.df a character type
target_meta.df$Dilution <- as.character(target_meta.df$Dilution)

#Make another target.df merged data frame for further use with tag-subtracted values and test samples only
target2.df <- merge(target_meta.df, norm_sub6.df, by.x = "Name", by.y ="row.names", all.y = TRUE, sort = FALSE)

#For seropositivity calculations, do once for Pf and once Pv, all antigens, all dilutions
Pf_antigens.df <- filter(target2.df, Plasmodium == "Pf")
Pf_antigens.df <- tibble::column_to_rownames(Pf_antigens.df, var="Name")
Pf_antigens.df <- Pf_antigens.df[,sapply(Pf_antigens.df, is.numeric)]

Pv_antigens.df <- filter(target2.df, Plasmodium == "Pv")
Pv_antigens.df <- tibble::column_to_rownames(Pv_antigens.df, var="Name")
Pv_antigens.df <- Pv_antigens.df[,sapply(Pv_antigens.df, is.numeric)]

#For the person_exposed calculation, only want the antigens at 1 dilution each.
#For the Sanger antigens Dilution = 0.5 
#For all other antigens, Dilution = 1
sub_Pf_antigens.df <- filter(target2.df, (Plasmodium == "Pf" & Dilution == "1") 
    | (Plasmodium == "Pf" & Source == "J. Rayner; WTSI" & Dilution == "0.5"))
sub_Pf_antigens.df <- tibble::column_to_rownames(sub_Pf_antigens.df, var="Name")
sub_Pf_antigens.df <- sub_Pf_antigens.df[,sapply(sub_Pf_antigens.df, is.numeric)]

sub_Pv_antigens.df <- filter(target2.df, (Plasmodium == "Pv" & Dilution == "1") 
    | (Plasmodium == "Pv" & Source == "J. Rayner; WTSI" & Dilution == "0.5"))
sub_Pv_antigens.df <- tibble::column_to_rownames(sub_Pv_antigens.df, var="Name")
sub_Pv_antigens.df <- sub_Pv_antigens.df[,sapply(sub_Pv_antigens.df, is.numeric)]

###Form a seropositivity matrix based on reactivity over a cutoff derived from sample buffer background.
#The threshold is the sample buffer mean + 3SD. Then take Log2 and normalize. 
sample_cutoff <- cor2_buffer_sample_mean + 3*cor2_buffer_sample_sd
log_sample_cutoff <- log2(sample_cutoff)
norm_sample_cutoff <- log_sample_cutoff - log_buffer_sample_mean

#Tailor the norm_sample_cutoff to remove excluded samples and control samples
buffer_cutoff.matrix <- as.matrix(norm_sample_cutoff)
rownames(buffer_cutoff.matrix, colnames(norm4.matrix))
sub_buf_cutoff.matrix <- as.matrix(buffer_cutoff.matrix[rownames(buffer_cutoff.matrix) %in% samples_test,])
sub_cutoff <- sub_buf_cutoff.matrix[(!rownames(sub_buf_cutoff.matrix) %in% samples_exclude),]

#Plot the sample cutoffs for samples included in analysis
png(filename = paste0(study, "_Buffer_Cutoffs.tif"), width = 5, height = 4, units = "in", res = 1200)
par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))
plot(sub_cutoff, pch='*', col = "blue", ylim=c(0,max(sub_cutoff)*1.25),
     ylab="Seropositivity Cutoff", xlab="Sample (Array)")

graphics.off()

#Then can apply the norm_sample_cutoff to each data frame of interest
#because the samples are not changing, only the antigens are changing. 

#Make seropositivity matrices for Pf and Pv separately 
SP_Pf.df <- t(apply(Pf_antigens.df, 1, function(x) ((x > sub_cutoff)+0)))
SP_Pv.df <- t(apply(Pv_antigens.df, 1, function(x) ((x > sub_cutoff)+0)))

sub_SP_Pf.df <- t(apply(sub_Pf_antigens.df, 1, function(x) ((x > sub_cutoff)+0)))
sub_SP_Pv.df <- t(apply(sub_Pv_antigens.df, 1, function(x) ((x > sub_cutoff)+0)))
  
###Create a threshold for overall target reactivity
#e.g. To be included in heatmaps and other analyses, perhaps targets should be reacted to by at least 5% of people?

#All Pf antigens - 80 out of 192 Pf targets are reactive in at least 5% of people
Pf_target_breadth <- rowSums(SP_Pf.df, na.rm=TRUE)
Pf_target_reactive <- Pf_target_breadth > (ncol(SP_Pf.df)/100)*5
cat(sum(Pf_target_reactive), "out of", nrow(SP_Pf.df), "Pf targets are reactive in at least 5% of people")

#All Pv antigens - 6 out of 63 Pv targets are reactive in at least 5% of people
Pv_target_breadth <- rowSums(SP_Pv.df, na.rm=TRUE)
Pv_target_reactive <- Pv_target_breadth > (ncol(SP_Pv.df)/100)*5
cat(sum(Pv_target_reactive), "out of", nrow(SP_Pv.df), "Pv targets are reactive in at least 5% of people")

###Create a threshold for overall person reactivity
#Similarly, perhaps unreactive individuals should be disincluded? Either way - this is informative.

#For person reactivity, we do not want to count multiple dilutions for each antigen 

#Sub Pf antigens - 1124 out of 1325 samples are reactive to at least 5% of Pf targets
Pf_person_breadth <- colSums(sub_SP_Pf.df, na.rm=TRUE)
Pf_person_exposed <- Pf_person_breadth > (nrow(sub_SP_Pf.df)/100)*5
cat(sum(Pf_person_exposed), "out of", ncol(sub_SP_Pf.df), "samples are reactive to at least 5% of Pf targets")

#Sub Pv antigens - 544 out of 1325 samples are reactive to at least 5% of Pv targets
Pv_person_breadth <- colSums(sub_SP_Pv.df, na.rm=TRUE)
Pv_person_exposed <- Pv_person_breadth > (nrow(sub_SP_Pv.df)/100)*5
cat(sum(Pv_person_exposed), "out of", ncol(sub_SP_Pv.df), "samples are reactive to at least 5% of Pv targets")

### Export matrix of data for reactive protein targets only (cutoff mean+3SD method)

# Includes control AND test samples
Pf.reactive.targets.matrix <- as.matrix(norm_sub4.df[Pf_target_reactive==TRUE,])
write.csv(Pf.reactive.targets.matrix, paste0(study,"Pf_reactive_targets_data_WG.csv")) 

Pv.reactive.targets.matrix <- as.matrix(norm_sub4.df[Pv_target_reactive==TRUE])
write.csv(Pv.reactive.targets.matrix, paste0(study,"Pv_reactive_targets_data_WG.csv")) 

#Pf plot of normalized data for each antigen organized by highest median
#Only include the data if the person is seropositive and an exposed person
exposed_SP_Pf.df <- SP_Pf.df[Pf_target_reactive==TRUE, Pf_person_exposed==TRUE]
reactive_Pf.df <- Pf_antigens.df[Pf_target_reactive==TRUE,Pf_person_exposed==TRUE]

#Export this matrix which only includes exposed individuals
write.csv(reactive_Pf.df, paste0(study,"Pf_exposed_data_WG.csv"))

#Plot the number of seropositive people by antigen, highest to lowest for Pf
SPpeople <- as.matrix(sort(rowSums(exposed_SP_Pf.df), decreasing = TRUE))
SPpeople <- as.data.frame(SPpeople)
SPpeople <- cbind(Target = rownames(SPpeople), SPpeople)
SPpeople$Target <- as.factor(SPpeople$Target)
#explicitly set factor levels to the correct order
SPpeople$Target <- factor(SPpeople$Target, levels = SPpeople$Target[order(-SPpeople$V1)])

png(filename = paste0(study, "_Pf_num_people_WG.tif"), width = 8, height = 4.5, units = "in", res = 1200)
par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))

ggplot(SPpeople, aes(x = Target, y = V1)) + theme_bw() + geom_bar(stat="identity") + ylab("Number of Seropositive Individuals") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 6))

graphics.off()
   
#Make a new data frame where seropositive values will be the number and otherwise it will be NA
SP_Pf_data.df <- data.frame(matrix(NA, nrow = nrow(reactive_Pf.df), ncol = ncol(reactive_Pf.df)))
rownames(SP_Pf_data.df) <- rownames(reactive_Pf.df)
colnames(SP_Pf_data.df) <- colnames(reactive_Pf.df)

for(b in 1:ncol(reactive_Pf.df)){
  for(a in 1:nrow(reactive_Pf.df)){
    if(exposed_SP_Pf.df[[a,b]]==1){
    SP_Pf_data.df[[a,b]] <- reactive_Pf.df[[a,b]] 
    }
  }
}
remove(a,b)

#then melt this data.frame with Na.rm = TRUE to organize for ggplot2
melt.Pf <- melt(as.matrix(SP_Pf_data.df), na.rm = TRUE)
colnames(melt.Pf) <- c("Target", "Sample", "Normalized")

#negative control data for Pf reactive targets, tag-subtracted data
neg_samples <-c(grep("Neg", colnames(norm4.matrix)))
neg_data <- norm_sub4.df[-rmsamp_all, neg_samples]
Pf_neg_data <- neg_data[Pf_target_reactive==TRUE,]
Pf_neg_mean <- as.matrix(rowMeans(Pf_neg_data))

#Add negative control data to the plot?

#Violin and Box Plots of data for reactive Pf antigens, sorted by highest median to lowest
ggplot(melt.Pf, aes(x=reorder(Target, -Normalized, FUN=median), y=Normalized)) + geom_violin()

png(filename = paste0(study, "_Pf_violin_WG.tif"), width = 5, height = 8, units = "in", res = 1200)
par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))

ggplot(melt.Pf, aes(x=reorder(Target, Normalized, FUN=median), y=Normalized)) + geom_violin() + coord_flip() + xlab("Target") + ylab("Normalized Log2(MFI)") + theme(text = element_text(size=9))

graphics.off()

png(filename = paste0(study, "_Pf_boxplot_WG.tif"), width = 5, height = 8, units = "in", res = 1200)
par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))

ggplot(melt.Pf, aes(x=reorder(Target, Normalized, FUN=median), y=Normalized)) + geom_boxplot(outlier.size = 0.3) + coord_flip() + xlab("Target") + ylab("Normalized Log2(MFI)") + theme(text = element_text(size=9))

graphics.off()

##### Repeat everything for Pf without the wheat germ antigens ######

#For seropositivity calculations, do once for Pf and once Pv, all antigens, all dilutions
Pf_no_WG.df <- filter(target2.df, Plasmodium == "Pf", !Expression_Tag=="Wheat Germ")
Pf_no_WG.df <- tibble::column_to_rownames(Pf_no_WG.df, var="Name")
Pf_no_WG.df <- Pf_no_WG.df[,sapply(Pf_no_WG.df, is.numeric)]

#For the person_exposed calculation, only want the antigens at 1 dilution each.
#For the Sanger antigens Dilution = 0.5 
#For all other antigens, Dilution = 1
sub_Pf_no_WG.df <- filter(target2.df, (Plasmodium == "Pf" & Dilution == "1" & !Expression_Tag=="Wheat Germ") 
                             | (Plasmodium == "Pf" & Source == "J. Rayner; WTSI" & Dilution == "0.5"))
sub_Pf_no_WG.df <- tibble::column_to_rownames(sub_Pf_no_WG.df, var="Name")
sub_Pf_no_WG.df <- sub_Pf_no_WG.df[,sapply(sub_Pf_no_WG.df, is.numeric)]

#Make seropositivity matrices for Pf , no wheat germ
SP_Pf_no_WG.df <- t(apply(Pf_no_WG.df, 1, function(x) ((x > sub_cutoff)+0)))
sub_SP_Pf_no_WG.df <- t(apply(sub_Pf_no_WG.df, 1, function(x) ((x > sub_cutoff)+0)))

###Create a threshold for overall target reactivity
#e.g. To be included in heatmaps and other analyses, perhaps targets should be reacted to by at least 5% of people?

#All Pf antigens - 73 out of 185 Pf targets (no WG) are reactive in at least 5% of people
Pf_target_breadth2 <- rowSums(SP_Pf_no_WG.df, na.rm=TRUE)
Pf_target_reactive2 <- Pf_target_breadth2 > (ncol(SP_Pf_no_WG.df)/100)*5
cat(sum(Pf_target_reactive2), "out of", nrow(SP_Pf_no_WG.df), "Pf targets (no WG) are reactive in at least 5% of people")

###Create a threshold for overall person reactivity
#Similarly, perhaps unreactive individuals should be disincluded? Either way - this is informative.

#For person reactivity, we do not want to count multiple dilutions for each antigen 

#Sub Pf antigens - 995 out of 1325 samples are reactive to at least 5% of Pf targets
Pf_person_breadth2 <- colSums(sub_SP_Pf_no_WG.df, na.rm=TRUE)
Pf_person_exposed2 <- Pf_person_breadth2 > (nrow(sub_SP_Pf_no_WG.df)/100)*5
cat(sum(Pf_person_exposed2), "out of", ncol(sub_SP_Pf_no_WG.df), "samples are reactive to at least 5% of Pf targets")

### Export matrix of data for reactive protein targets only (cutoff mean+3SD method)

# Includes control AND test samples, no wheat germ AG
Pf.reactive.targets.matrix <- as.matrix(norm_sub4.df[Pf_target_reactive2==TRUE,])
write.csv(Pf.reactive.targets.matrix, paste0(study,"Pf_reactive_targets_data.csv")) 

###below this line, stopped making separate names for variables for w/o WG antigens.

#Pf plot of normalized data for each antigen organized by highest median
#Only include the data if the person is seropositive and an exposed person
exposed_SP_Pf.df <- SP_Pf_no_WG.df[Pf_target_reactive2==TRUE, Pf_person_exposed2==TRUE]
reactive_Pf.df <- Pf_no_WG.df[Pf_target_reactive2==TRUE,Pf_person_exposed2==TRUE]

#Export this matrix which only includes exposed individuals
write.csv(reactive_Pf.df, paste0(study,"Pf_exposed_data.csv"))

#Plot the number of seropositive people by antigen, highest to lowest for Pf
SPpeople <- as.matrix(sort(rowSums(exposed_SP_Pf.df), decreasing = TRUE))
SPpeople <- as.data.frame(SPpeople)
SPpeople <- cbind(Target = rownames(SPpeople), SPpeople)
SPpeople$Target <- as.factor(SPpeople$Target)
#explicitly set factor levels to the correct order
SPpeople$Target <- factor(SPpeople$Target, levels = SPpeople$Target[order(-SPpeople$V1)])

png(filename = paste0(study, "_Pf_num_people.tif"), width = 8, height = 4.5, units = "in", res = 1200)
par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))

ggplot(SPpeople, aes(x = Target, y = V1)) + theme_bw() + geom_bar(stat="identity") + ylab("Number of Seropositive Individuals") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 6))

graphics.off()

#Make a new data frame where seropositive values will be the number and otherwise it will be NA
SP_Pf_data.df <- data.frame(matrix(NA, nrow = nrow(reactive_Pf.df), ncol = ncol(reactive_Pf.df)))
rownames(SP_Pf_data.df) <- rownames(reactive_Pf.df)
colnames(SP_Pf_data.df) <- colnames(reactive_Pf.df)

for(b in 1:ncol(reactive_Pf.df)){
  for(a in 1:nrow(reactive_Pf.df)){
    if(exposed_SP_Pf.df[[a,b]]==1){
      SP_Pf_data.df[[a,b]] <- reactive_Pf.df[[a,b]] 
    }
  }
}
remove(a,b)

#then melt this data.frame with Na.rm = TRUE to organize for ggplot2
melt.Pf <- melt(as.matrix(SP_Pf_data.df), na.rm = TRUE)
colnames(melt.Pf) <- c("Target", "Sample", "Normalized")

#negative control data for Pf reactive targets, tag-subtracted data
neg_samples <-c(grep("Neg", colnames(norm4.matrix)))
neg_data <- norm_sub4.df[-rmsamp_all, neg_samples]
Pf_neg_data <- neg_data[Pf_target_reactive2==TRUE,]
Pf_neg_mean <- as.matrix(rowMeans(Pf_neg_data))

#Add negative control data to the plot?

#Violin and Box Plots of data for reactive Pf antigens, sorted by highest median to lowest
ggplot(melt.Pf, aes(x=reorder(Target, -Normalized, FUN=median), y=Normalized)) + geom_violin()

png(filename = paste0(study, "_Pf_violin.tif"), width = 5, height = 8, units = "in", res = 1200)
par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))

ggplot(melt.Pf, aes(x=reorder(Target, Normalized, FUN=median), y=Normalized)) + geom_violin() + coord_flip() + xlab("Target") + ylab("Normalized Log2(MFI)") + theme(text = element_text(size=9))

graphics.off()

png(filename = paste0(study, "_Pf_boxplot.tif"), width = 5, height = 8, units = "in", res = 1200)
par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))

ggplot(melt.Pf, aes(x=reorder(Target, Normalized, FUN=median), y=Normalized)) + geom_boxplot(outlier.size = 0.3) + coord_flip() + xlab("Target") + ylab("Normalized Log2(MFI)") + theme(text = element_text(size=9))

graphics.off()

######### Repeating everything for Pf again below for malarial centre abstract!!!!!!
#### not changing names of data frames for this
# do not export the data

##### Repeat everything for Pf without Particular Sources and  without the wheat germ antigens ######

#because we do want MSP1-19 from AA Holder, change the source column for that antigen only
target2.df$Source[target2.df$Name == "PfMSP1_19"] <- "A.A. Holder MSP1-19" 

#For seropositivity calculations, do once for Pf and once Pv, all antigens, all dilutions
Pf_no_WG.df <- filter(target2.df, Plasmodium == "Pf", !Expression_Tag=="Wheat Germ", 
    !Source == "J. Beeson plus; Burnet", !Source == "A.A. Holder; Francis Crick", !Source =="C. Kocken; BPRC",
    !Source == "A. Jensen; CMP", !Source == "R. Coppel; Monash", !Source == "S. Draper; Jenner", !Source == "B. Arca; La Sapienza")

#only want the antigens at 1 dilution each.
#For the Sanger antigens Dilution = 0.5 
#For all other antigens, Dilution = 1
sub_Pf_no_WG.df <- filter(Pf_no_WG.df, Dilution == "1" | (Source == "J. Rayner; WTSI" & Dilution == "0.5"))
sub_Pf_no_WG.df <- tibble::column_to_rownames(sub_Pf_no_WG.df, var="Name")
sub_Pf_no_WG.df <- sub_Pf_no_WG.df[,sapply(sub_Pf_no_WG.df, is.numeric)]

#Make seropositivity matrices for Pf , no wheat germ
sub_SP_Pf_no_WG.df <- t(apply(sub_Pf_no_WG.df, 1, function(x) ((x > sub_cutoff)+0)))

###Create a threshold for overall target reactivity
#e.g. To be included in heatmaps and other analyses, perhaps targets should be reacted to by at least 5% of people?

#Selected Pf antigens - 57 out of 111 Pf targets (no WG) are reactive in at least 5% of people
Pf_target_breadth2 <- rowSums(sub_Pf_no_WG.df, na.rm=TRUE)
Pf_target_reactive2 <- Pf_target_breadth2 > (ncol(sub_Pf_no_WG.df)/100)*5
cat(sum(Pf_target_reactive2), "out of", nrow(sub_Pf_no_WG.df), "Pf targets (no WG) are reactive in at least 5% of people")

###Create a threshold for overall person reactivity
#Similarly, perhaps unreactive individuals should be disincluded? Either way - this is informative.

#For person reactivity, we do not want to count multiple dilutions for each antigen 

#Sub Pf antigens - selected targets - 979 out of 1325 samples are reactive to at least 5% of Pf targets
Pf_person_breadth2 <- colSums(sub_SP_Pf_no_WG.df, na.rm=TRUE)
Pf_person_exposed2 <- Pf_person_breadth2 > (nrow(sub_SP_Pf_no_WG.df)/100)*5
cat(sum(Pf_person_exposed2), "out of", ncol(sub_SP_Pf_no_WG.df), "samples are reactive to at least 5% of Pf targets")

#Pf plot of normalized data for each antigen organized by highest median
#Only include the data if the person is seropositive and an exposed person
exposed_SP_Pf.df <- sub_SP_Pf_no_WG.df[Pf_target_reactive2==TRUE, Pf_person_exposed2==TRUE]
reactive_Pf.df <- sub_Pf_no_WG.df[Pf_target_reactive2==TRUE,Pf_person_exposed2==TRUE]

#Plot the number of seropositive people by antigen, highest to lowest for Pf
SPpeople <- as.matrix(sort(rowSums(exposed_SP_Pf.df), decreasing = TRUE))
SPpeople <- as.data.frame(SPpeople)
SPpeople <- cbind(Target = rownames(SPpeople), SPpeople)
SPpeople$Target <- as.factor(SPpeople$Target)
#explicitly set factor levels to the correct order
SPpeople$Target <- factor(SPpeople$Target, levels = SPpeople$Target[order(-SPpeople$V1)])

png(filename = paste0(study, "_Pf_num_people_MCabs.tif"), width = 8, height = 4.5, units = "in", res = 1200)
par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))

ggplot(SPpeople, aes(x = Target, y = V1)) + theme_bw() + geom_bar(stat="identity") + ylab("Number of Seropositive Individuals") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 6))

graphics.off()

#Make a new data frame where seropositive values will be the number and otherwise it will be NA
SP_Pf_data.df <- data.frame(matrix(NA, nrow = nrow(reactive_Pf.df), ncol = ncol(reactive_Pf.df)))
rownames(SP_Pf_data.df) <- rownames(reactive_Pf.df)
colnames(SP_Pf_data.df) <- colnames(reactive_Pf.df)

for(b in 1:ncol(reactive_Pf.df)){
  for(a in 1:nrow(reactive_Pf.df)){
    if(exposed_SP_Pf.df[[a,b]]==1){
      SP_Pf_data.df[[a,b]] <- reactive_Pf.df[[a,b]] 
    }
  }
}
remove(a,b)

#then melt this data.frame with Na.rm = TRUE to organize for ggplot2
melt.Pf <- melt(as.matrix(SP_Pf_data.df), na.rm = TRUE)
colnames(melt.Pf) <- c("Target", "Sample", "Normalized")

#negative control data for Pf reactive targets, tag-subtracted data
neg_samples <-c(grep("Neg", colnames(norm4.matrix)))
neg_data <- norm_sub4.df[-rmsamp_all, neg_samples]
Pf_neg_data <- neg_data[Pf_target_reactive2==TRUE,]
Pf_neg_mean <- as.matrix(rowMeans(Pf_neg_data))

#Add negative control data to the plot?

#Violin and Box Plots of data for reactive Pf antigens, sorted by highest median to lowest
ggplot(melt.Pf, aes(x=reorder(Target, -Normalized, FUN=median), y=Normalized)) + geom_violin()

png(filename = paste0(study, "_Pf_violin_MCabs.tif"), width = 5, height = 8, units = "in", res = 1200)
par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))

ggplot(melt.Pf, aes(x=reorder(Target, Normalized, FUN=median), y=Normalized)) + geom_violin() + coord_flip() + xlab("Target") + ylab("Normalized Log2(MFI)") + theme(text = element_text(size=9))

graphics.off()

png(filename = paste0(study, "_Pf_boxplot_MCabs.tif"), width = 5, height = 8, units = "in", res = 1200)
par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))

ggplot(melt.Pf, aes(x=reorder(Target, Normalized, FUN=median), y=Normalized)) + geom_boxplot(outlier.size = 0.3) + coord_flip() + xlab("Target") + ylab("Normalized Log2(MFI)") + theme(text = element_text(size=9))

graphics.off()



#######Same thing for Pv Antigens########

#Only include the data if the person is seropositive and it's an exposed individual
exposed_SP_Pv.df <- SP_Pv.df[Pv_target_reactive==TRUE, Pv_person_exposed==TRUE]
reactive_Pv.df <- Pv_antigens.df[Pv_target_reactive==TRUE,Pv_person_exposed==TRUE]

#Export this matrix which only includes exposed individuals
write.csv(reactive_Pv.df, paste0(study,"Pv_exposed_data.csv"))

#Plot the number of seropositive people by antigen, highest to lowest for Pf
SPpeoplePv <- as.matrix(sort(rowSums(exposed_SP_Pv.df), decreasing = TRUE))
SPpeoplePv <- as.data.frame(SPpeoplePv)
SPpeoplePv <- cbind(Target = rownames(SPpeoplePv), SPpeoplePv)
SPpeoplePv$Target <- as.factor(SPpeoplePv$Target)
#explicitly set factor levels to the correct order
SPpeoplePv$Target <- factor(SPpeoplePv$Target, levels = SPpeoplePv$Target[order(-SPpeoplePv$V1)])

png(filename = paste0(study, "_Pv_num_people.tif"), width = 3, height = 4.5, units = "in", res = 1200)
par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))

ggplot(SPpeoplePv, aes(x = Target, y = V1)) + theme_bw() + geom_bar(stat="identity") + ylab("Number of Seropositive Individuals") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 6))

graphics.off()

#Make a new data frame where seropositive values will be the number and otherwise it will be NA
SP_Pv_data.df <- data.frame(matrix(NA, nrow = nrow(reactive_Pv.df), ncol = ncol(reactive_Pv.df)))
rownames(SP_Pv_data.df) <- rownames(reactive_Pv.df)
colnames(SP_Pv_data.df) <- colnames(reactive_Pv.df)

for(b in 1:ncol(reactive_Pv.df)){
  for(a in 1:nrow(reactive_Pv.df)){
    if(exposed_SP_Pv.df[[a,b]]==1){
      SP_Pv_data.df[[a,b]] <- reactive_Pv.df[[a,b]] 
    }
  }
}
remove(a,b)

#negative control data for Pv reactive targets, tag-subtracted data
Pv_neg_data <- neg_data[Pv_target_reactive==TRUE,]
Pv_neg_mean <- as.matrix(rowMeans(Pv_neg_data))

#then melt this data.frame with Na.rm = TRUE to organize for ggplot2
melt.Pv <- melt(as.matrix(SP_Pv_data.df), na.rm = TRUE)
colnames(melt.Pv) <- c("Target", "Sample", "Normalized")

#Violin and Box Plots of data for reactive Pv antigens, sorted by highest median to lowest
ggplot(melt.Pv, aes(x=reorder(Target, -Normalized, FUN=median), y=Normalized)) + geom_violin()

png(filename = paste0(study, "_Pv_violin.tif"), width = 5, height = 4, units = "in", res = 1200)

ggplot(melt.Pv, aes(x=reorder(Target, -Normalized, FUN=median), y=Normalized)) + geom_violin() + 
  xlab("Target") + ylab("Normalized Log2(MFI)") + 
  theme(text = element_text(size=10), axis.text.x = element_text(size = 8, angle = 45, hjust = 1, vjust = 1))

graphics.off()

png(filename = paste0(study, "_Pv_boxplot.tif"), width = 5, height = 4, units = "in", res = 1200)

ggplot(melt.Pv, aes(x=reorder(Target, -Normalized, FUN=median), y=Normalized)) + geom_boxplot(outlier.size = 0.3) + 
  xlab("Target") + ylab("Normalized Log2(MFI)") + 
  theme(text = element_text(size=10), axis.text.x = element_text(size = 8, angle = 45, hjust = 1, vjust = 1))

graphics.off()

save.image("SangerAnalysis.RData")

#Old code - not using now
### Plot of geometric mean vs target, ranked from highest to lowest

# # Only using data from reactive targets and reactive test samples
# reactive.matrix <- reactive.targets.matrix[,person_exposed]
# 
# #Calculate geometric mean and geometric SD for each antigen
# #I have not done anything with the seropositivity / seronegativity here
# mean_targets <- rowMeans(reactive.matrix)
# sd_targets <- apply(reactive.matrix, 1, sd)
# 
# target_data <- data.frame(mean_targets, sd_targets, neg_mean)
# target_data <- target_data[order(-mean_targets),]
# 
# sorted_mean <- mean_targets[order(-mean_targets)]
# 
# 
# #barPlot - this still looks terrible 
# png(filename = paste0(study, "_GeoMean_Reactive_Targets.tif"), width = 5, height = 4, units = "in", res = 1200)
# par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))
# 
# barplot(sorted_mean, pch='*', col = "blue", ylim=c(0,max(mean_targets)*1.25))
# 
# axis(2, ylab="Geometric Mean Log2(Target MFI/Buffer MFI)")
# axis(1, cex.axis = 0.4, labels = as.character(row.names(target_data)), at = 1:length(sorted_mean), xlab="Antigen")
# 
# graphics.off()
# 
# ### Plot of number of seropositive individuals for each sanger antigen 
# reactive_seroposSD.matrix <- seroposSD.matrix[target_reactive==TRUE, person_exposed]
# rownames(reactive_seroposSD.matrix) <- rownames(reactive.targets.matrix)
# sub_sanger_antigens <- c(grep("(s)", rownames(reactive_seroposSD.matrix), fixed = TRUE))
# sanger_seroposSD.matrix <- reactive_seroposSD.matrix[sub_sanger_antigens,]
# 
# sanger_seroposSD.df <- tibble::rownames_to_column(sanger_seroposSD.matrix)
# 
# #Sum of people positive for each reactive sanger antigen, sorted highest to lowest
# sanger_sums <- sort(rowSums(sanger_seroposSD.matrix), decreasing = TRUE)
# #check plot in R
# barplot(sanger_sums)
# 
# #This plot still doesn't look good, I just didn't finish and make it quickly in excel to move on
# png(filename = paste0(study,"_targets_seropos_sums.tif"), width = 5.5, height = 4, units = "in", res = 600)
# par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1), bty = "o", 
#     mgp = c(2, 0.5, 0), cex.main = 0.7, cex.axis = 0.7, cex.lab = 1, xpd=NA, las=1)
# 
# barplot(sanger_sums, xlab="Target", add=FALSE, col = "darkblue", ylim=c(0,max(sanger_sums)*1.25))
# title(main = "Number of People Seropositive for each Reactive Antigen", adj=0)
# title(ylab="Number of People", line=2.7)
# 
# graphics.off()
# 
# qplot(t(sanger_seroposSD.matrix), geom="histogram")
# 
# 
# 
# 
# 

