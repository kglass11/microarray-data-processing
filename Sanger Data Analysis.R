#Sanger Data Analysis Script
#Katie Glass

#####################################
############DATA ANALYSIS############
#####################################

#If you haven't just continued from the processing script, run:
#load(file="SangerAfterProcessing.RData")

#install.packages("reshape2")

library(reshape2)

###Seropositivity and Reactivity Thresholds###

#For this section, using normalized data with negative values set to 0 (norm4.matrix)

#Before doing any further analysis, we have to get rid of samples or targets that we are no longer 
#interested in.
#E.g. If control individuals are in our analysis, they will affect mixture model based cut-offs
#E.g. If control targets are still in our analysis, they will muck up our protein breadth estimates
#KG - for now this still includes the same protein target at different dilutions
#This means we have to subset the data, so some earlier annotations will from here on be wrong 
#(e.g. index_sample will no longer equal 96)

#Assign sample type names, to identify control and test samples (logical)
samples_test <- sample_meta_f.df$sample_type=="test"
samples_control <- sample_meta_f.df$sample_type=="control"

#Define a list of targets to be removed from further analysis (controls)
rmsamp_all <- unique(c(targets_blank, targets_buffer, targets_ref, targets_std, high_targets_disinclude))

#Remove samples that should be excluded
norm_sub.matrix <- norm4.matrix[,(!colnames(norm4.matrix) %in% samples_exclude)]

#Remove control protein targets and control samples for seropositivity calculations
norm_sub2.matrix <- norm_sub.matrix[-rmsamp_all, samples_test]
samples_sub.df <- sample_meta_f.df[samples_test,]

#Replace current target names with original target names now that control targets are removed
#might be useful to merge this instead with the target dataframe?
norm_sub3.df <- merge(norm_sub2.matrix, annotation_targets.df, by ="row.names", sort = FALSE)
norm_sub3.df <- tibble::column_to_rownames(norm_sub3.df, var="Row.names")
row.names(norm_sub3.df) <- norm_sub3.df$Name
norm_sub4.df <- norm_sub3.df[,1:ncol(norm_sub2.matrix)]

#Merge with target metadata to filter based on expression tag etc.
target.df <- merge(target_meta.df, norm_sub4.df, by.x = "Name", by.y ="row.names", all.y = TRUE, sort = FALSE)
  
### Subtracting Protein Tag Signal from tagged antigens - only for norm4.matrix (no negative values)
#Prepare data frame with GST tagged proteins only for subtraction
GST_antigens.df <- filter(target.df, Expression_Tag == "GST" | Expression_Tag == "GST/His")
GST_antigens.df <- tibble::column_to_rownames(GST_antigens.df, var="Name")
GST_antigens.df <- GST_antigens.df[,sapply(GST_antigens.df, is.numeric)]

#Subtract GST signal from GST tagged proteins
GST <- c(grep("GST", rownames(norm_sub4.df), fixed = TRUE))

sub_GST_antigens.df <- data.frame(matrix(0, nrow = nrow(GST_antigens.df), ncol = ncol(GST_antigens.df)))
rownames(sub_GST_antigens.df) <- rownames(GST_antigens.df)
colnames(sub_GST_antigens.df) <- colnames(GST_antigens.df)

#All the NAs have been removed already by excluding targets and samples
for(b in 1:ncol(GST_antigens.df)){
  for(a in 1:nrow(GST_antigens.df)){
  # when the GST value is positive only, subtract GST value, 
  #otherwise want to leave as whatever the value was before (because GST was at or below buffer mean for that sample)
    if (norm_sub4.df[GST,b] > 0){
        #calculate difference in original MFI form (not log2)
        sub_GST_antigens.df[a,b] <- 2^GST_antigens.df[a,b] - 2^norm_sub4.df[GST,b]
        #can only do log2 if the difference is greater than 0, otherwise set to 0 (normalized log2 value is not above buffer)
        if (sub_GST_antigens.df[a,b] > 0) {
          sub_GST_antigens.df[a,b] <- log2(sub_GST_antigens.df[a,b])
          #if the log2 tag-subtracted value is negative, means normalized value is below buffer mean,
          #so also need to set those negatives to 0 again.
            if(sub_GST_antigens.df[a,b] < 0){
              sub_GST_antigens.df[a,b] <- 0
            }
        } else { 
          sub_GST_antigens.df[a,b] <- 0
          }
    } else {
      sub_GST_antigens.df[a,b] <- GST_antigens.df[a,b]
      }
  }
}
remove(a,b)

#Prepare CD4-tagged only data frame (Sanger antigens)
CD4_antigens.df <- filter(target.df, Expression_Tag == "CD4")
CD4_antigens.df <- tibble::column_to_rownames(CD4_antigens.df, var="Name")
CD4_antigens.df <- CD4_antigens.df[,sapply(CD4_antigens.df, is.numeric)]

#Subtract CD4 from CD4 tagged antigens
CD4 <- c(grep("CD4", rownames(norm_sub4.df), fixed = TRUE))

sub_CD4_antigens.df <- data.frame(matrix(0, nrow = nrow(CD4_antigens.df), ncol = ncol(CD4_antigens.df)))
rownames(sub_CD4_antigens.df) <- rownames(CD4_antigens.df)
colnames(sub_CD4_antigens.df) <- colnames(CD4_antigens.df)

#All the NAs have been removed already by excluding targets and samples
for(b in 1:ncol(CD4_antigens.df)){
  for(a in 1:nrow(CD4_antigens.df)){
    # when the CD4 value is positive only, subtract CD4 value, 
    #otherwise want to leave as whatever the value was before (because CD4 was at or below buffer mean for that sample)
    if (norm_sub4.df[CD4,b] > 0){
      #calculate difference in original MFI form (not log2)
      sub_CD4_antigens.df[a,b] <- 2^CD4_antigens.df[a,b] - 2^norm_sub4.df[CD4,b]
      #can only do log2 if the difference is greater than 0, otherwise set to 0 (normalized log2 value is not above buffer)
      if (sub_CD4_antigens.df[a,b] > 0) {
        sub_CD4_antigens.df[a,b] <- log2(sub_CD4_antigens.df[a,b])
        #if the log2 tag-subtracted value is negative, means normalized value is below buffer mean,
        #so also need to set those negatives to 0 again.
        if(sub_CD4_antigens.df[a,b] < 0){
          sub_CD4_antigens.df[a,b] <- 0
        }
      } else { 
        sub_CD4_antigens.df[a,b] <- 0
      }
    } else {
      sub_CD4_antigens.df[a,b] <- CD4_antigens.df[a,b]
    }
  }
}
remove(a,b)

### Plot CD4 and GST signals
GST_val <- c(as.matrix(norm_sub4.df[GST,]))
CD4_val <- c(as.matrix(norm_sub4.df[CD4,]))

png(filename = paste0(study, "_Tags_GST_CD4.tif"), width = 5, height = 7.5, units = "in", res = 1200)
par(mfrow=c(2,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))
plot(GST_val, pch='*', col = "blue", ylim=c(0,max(GST_val)*1.25), main = "GST",
     ylab="Normalized log2(MFI)", xlab="Sample (Array)")

plot(CD4_val, pch='*', col = "darkblue", ylim=c(0,max(CD4_val)*1.25), main = "CD4",
     ylab="Normalized log2(MFI)", xlab="Sample (Array)")

#print text on the plots for number of samples where GST and CD4 are above buffer
mtext(paste("Total Samples with GST > 0:", round(sum(norm_sub4.df[GST,] > 0), digits=2), 
  "(", round(sum(norm_sub4.df[GST,] > 0)/length(GST_val)*100, digits=2), "%)"), side=1, cex=0.8, line=0.5, outer=TRUE, xpd=NA, adj=0)
mtext(paste("Total Samples with CD4 > 0:", round(sum(norm_sub4.df[CD4,] > 0), digits=2), 
  "(", round(sum(norm_sub4.df[CD4,] > 0)/length(CD4_val)*100, digits=2), "%)"), side=1, cex=0.8, line=1.5, outer=TRUE, xpd=NA, adj=0)

graphics.off()

#Make another data frame where the tagged protein values are replaced by their subtracted values
#filter out the GST and CD4 proteins
no_tags.df <- filter(target.df, !(Expression_Tag == "CD4" | Expression_Tag == "GST" | Expression_Tag == "GST/His"))
no_tags.df <- tibble::column_to_rownames(no_tags.df, var="Name")
no_tags.df <- no_tags.df[,sapply(no_tags.df, is.numeric)]
#then rbind the GST and the CD4 data frames to that one. The order of the targets
#shouldn't matter anymore
norm_sub5.df <- rbind(no_tags.df, sub_GST_antigens.df, sub_CD4_antigens.df)

#Make another target.df merged data frame for further use with tag-subtracted values
target2.df <- merge(target_meta.df, norm_sub5.df, by.x = "Name", by.y ="row.names", all.y = TRUE, sort = FALSE)

#For seropositivity calculations, do once for Pf and once Pv, all antigens, regardless of dilution
#In future, can add a column to target metadata for "Dilution" and sort based on Dilution == "1"
Pf_antigens.df <- filter(target2.df, Plasmodium == "Pf")
Pf_antigens.df <- tibble::column_to_rownames(Pf_antigens.df, var="Name")
Pf_antigens.df <- Pf_antigens.df[,sapply(Pf_antigens.df, is.numeric)]

Pv_antigens.df <- filter(target2.df, Plasmodium == "Pv")
Pv_antigens.df <- tibble::column_to_rownames(Pv_antigens.df, var="Name")
Pv_antigens.df <- Pv_antigens.df[,sapply(Pv_antigens.df, is.numeric)]

###Form a seropositivity matrix based on reactivity over a cutoff derived from sample buffer background.
#The threshold is the sample buffer mean + 3SD. 
sample_cutoff <- cor2_buffer_sample_mean + 3*cor2_buffer_sample_sd
log_sample_cutoff <- log2(sample_cutoff)
norm_sample_cutoff <- log_sample_cutoff - log_buffer_sample_mean

#This seropositivity matrix includes ALL non-control antigens and all non-excluded "test" samples
#*** This is still using the data before switching to subtracted antigens! need to change strategy on this part
seroposSD_temp.matrix <- t(apply(norm4.matrix, 1, function(x) (x > norm_sample_cutoff)+0))
seroposSD_temp.matrix <- seroposSD_temp.matrix[,(!colnames(seroposSD_temp.matrix) %in% samples_exclude)]
seroposSD.matrix <- seroposSD_temp.matrix[-rmsamp_all, samples_test]

#Make seropositivity matrices for Pf and Pv separately 

# *** We might switch to using neg populations run on microarray slides for the cutoffs


###Create a threshold for overall target and person reactivity
#e.g. To be included in heatmaps and other analyses, perhaps targets should be reacted to by at least 5% of people?
#Similarly, perhaps unreactive individuals should be disincluded? Either way - this is informative.
target_breadth <- rowSums(seroposSD.matrix, na.rm=TRUE)
target_reactive <- target_breadth > (ncol(seroposSD.matrix)/100)*5
cat(sum(target_reactive), "out of", nrow(seroposSD.matrix), "protein targets are reactive in at least 5% of people")

person_breadth <- colSums(seroposSD.matrix, na.rm=TRUE)
person_exposed <- person_breadth > (nrow(seroposSD.matrix)/100)*5
cat(sum(person_exposed), "out of", ncol(seroposSD.matrix), "samples are reactive to at least 5% of proteins")

### Export matrix of data for reactive protein targets only (cutoff mean+3SD method)

# Includes control and test samples but not excluded samples
reactive.targets.matrix <- as.matrix(norm_sub4.df[target_reactive==TRUE,])
write.csv(reactive.targets.matrix, paste0(study,"_reactive_targets_data.csv")) 

### Plot of geometric mean vs target, ranked from highest to lowest

# Only using data from reactive targets and reactive test samples
reactive.matrix <- reactive.targets.matrix[,person_exposed]

#negative control data for reactive targets
neg_samples <-c(grep("Neg", colnames(norm4.matrix)))
neg_data <- norm4.matrix[-rmsamp_all, neg_samples]
neg_data <- neg_data[target_reactive==TRUE,]
neg_mean <- rowMeans(neg_data)

#Calculate geometric mean and geometric SD for each antigen
#I have not done anything with the seropositivity / seronegativity here
mean_targets <- rowMeans(reactive.matrix)
sd_targets <- apply(reactive.matrix, 1, sd)

target_data <- data.frame(mean_targets, sd_targets, neg_mean)
target_data <- target_data[order(-mean_targets),]

sorted_mean <- mean_targets[order(-mean_targets)]

#Violin Plot for sorted means! :) 

#organize the data for ggplot2 - using reactive.matrix for now
melt.reactive <- melt(reactive.matrix)
colnames(melt.reactive) <- c("Target", "Sample", "Normalized")

#sort data by mean - this isn't working yet
sort.melt.reactive <- with(melt.reactive, reorder(Normalized, Target, FUN = mean))

ggplot(melt.reactive, aes(x=Target, y=Normalized)) + geom_violin(trim = FALSE)

#1st 11 antigens only - this is not the final data just a test to see the plots
melt.11 <- melt(as.matrix(reactive.matrix[1:11,]))
colnames(melt.11) <- c("Target", "Sample", "Normalized") 

png(filename = paste0(study, "_11_test.tif"), width = 8, height = 4, units = "in", res = 1200)
par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))

ggplot(melt.11, aes(x=Target, y=Normalized)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

graphics.off()

#barPlot - this still looks terrible 
png(filename = paste0(study, "_GeoMean_Reactive_Targets.tif"), width = 5, height = 4, units = "in", res = 1200)
par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))

barplot(sorted_mean, pch='*', col = "blue", ylim=c(0,max(mean_targets)*1.25))

axis(2, ylab="Geometric Mean Log2(Target MFI/Buffer MFI)")
axis(1, cex.axis = 0.4, labels = as.character(row.names(target_data)), at = 1:length(sorted_mean), xlab="Antigen")

graphics.off()

### Plot of number of seropositive individuals for each sanger antigen 
reactive_seroposSD.matrix <- seroposSD.matrix[target_reactive==TRUE, person_exposed]
rownames(reactive_seroposSD.matrix) <- rownames(reactive.targets.matrix)
sub_sanger_antigens <- c(grep("(s)", rownames(reactive_seroposSD.matrix), fixed = TRUE))
sanger_seroposSD.matrix <- reactive_seroposSD.matrix[sub_sanger_antigens,]

sanger_seroposSD.df <- tibble::rownames_to_column(sanger_seroposSD.matrix)

#Sum of people positive for each reactive sanger antigen, sorted highest to lowest
sanger_sums <- sort(rowSums(sanger_seroposSD.matrix), decreasing = TRUE)
#check plot in R
barplot(sanger_sums)

#This plot still doesn't look good, I just didn't finish and make it quickly in excel to move on
png(filename = paste0(study,"_targets_seropos_sums.tif"), width = 5.5, height = 4, units = "in", res = 600)
par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1), bty = "o", 
    mgp = c(2, 0.5, 0), cex.main = 0.7, cex.axis = 0.7, cex.lab = 1, xpd=NA, las=1)

barplot(sanger_sums, xlab="Target", add=FALSE, col = "darkblue", ylim=c(0,max(sanger_sums)*1.25))
title(main = "Number of People Seropositive for each Reactive Antigen", adj=0)
title(ylab="Number of People", line=2.7)

graphics.off()

qplot(t(sanger_seroposSD.matrix), geom="histogram")






