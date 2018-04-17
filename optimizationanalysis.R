#optimization data analysis 
#KG

#"I:/Drakeley Group/Protein microarrays/Experiments/270717 Optimisation"
#
setwd("/Users/Katie/Desktop/R files from work/270717 Optimisation")
getwd()

# install.packages("lme4")
# install.packages("broom")
# install.packages("multcomp")
# install.packages("lsmeans")
# install.packages("emmeans")
# install.packages("pbkrtest")
# install.packages("phia")
# install.packages("lmerTest")
# install.packages("multcompView")
# install.packages("MuMIn")
# install.packages("NeatMap")

library(broom)
library(lme4)
library(multcomp)
library(emmeans)
library(pbkrtest)
library(phia)
library(lmerTest)
library(multcompView)
library(MuMIn)

require("gtools")

library(limma)
library(contrast)
library(beeswarm)
library(mixtools)
library(gplots)
library(ggplot2)
library(gcookbook)
library(dplyr)
library(reshape2)
library(NeatMap)

load("OptimizationLMMready.RData")

#this is the file from Tate's data processing
#load("ATT93447.RData")

#this is with older method of patrick's rule and original LMM
#load("OptimizationAfterPairwise.RData")

#import target metadata file
target_file <- "Opto Target Metadata.csv" 
target_meta.df <- read.csv(target_file, header=T, na.strings = " ", check.names = FALSE, stringsAsFactors = FALSE)

# set negative normalized values to zero - replicates are stored in trans.norm.rep1 and trans.norm.rep2
rep1.matrix <- as.matrix(trans.norm.rep1)

for (i in 1:length(rep1.matrix))
{
  if (is.na(rep1.matrix[[i]]) | rep1.matrix[[i]] > 0) 
  {
    i = i+1
  } else if(rep1.matrix[[i]] < 0) 
  {
    rep1.matrix[[i]] <- 0
  }
}

rep2.matrix <- as.matrix(trans.norm.rep2)

for (i in 1:length(rep2.matrix))
{
  if (is.na(rep2.matrix[[i]]) | rep2.matrix[[i]] > 0) 
  {
    i = i+1
  } else if(rep2.matrix[[i]] < 0) 
  {
    rep2.matrix[[i]] <- 0
  }
}

### Average duplicates
  trans.norm.avg <- matrix(nrow = nrow(trans.norm.rep1), ncol = ncol(trans.norm.rep1))
  colnames(trans.norm.avg) = colnames(trans.norm.rep1)
  rownames(trans.norm.avg) = rownames(trans.norm.rep1)
  
  trans.norm.avg <- log2((2^rep1.matrix + 2^rep2.matrix)/2)
  
above2 <- which(trans.norm.avg > 2)

### Check for deviant technical replicates, automatically exclude (set to NA)
# Use Patrick's formula for ELISA to compare replicates within one array
# if rep1 or rep2 is more than 1.5 times rep2 or rep1, respectively, exclude that pair
# Added: if either rep1 or rep2 is less than 0.2, then don't apply the rule
# for(k in 1:ncol(rep1.matrix)){
#    for(j in 1:nrow(rep1.matrix)){
#      if(is.na(rep1.matrix[j,k]) | is.na(rep2.matrix[j,k]) | rep1.matrix[j,k]<0.2 | rep2.matrix[j,k]<0.2){
#        j+1
#      } else if (rep1.matrix[j,k] > (log2(1.5) + rep2.matrix[j,k]) | (rep2.matrix[j,k] > (log2(1.5) + rep1.matrix[j,k])) == TRUE) 
#      {
#        trans.norm.avg[j,k] <- NA
#      }
#    }
#  }
# remove(j,k)

#adjusting Patrick's rule for this data: 
#check how many NA values there are with above method
# sum(is.na(trans.norm.avg))/length(c(trans.norm.avg))
#this is actually only 5.54% NA haha

#try changing the rule to only apply to those less than 2
# for(k in 1:ncol(rep1.matrix)){
#   for(j in 1:nrow(rep1.matrix)){
#     if(is.na(rep1.matrix[j,k]) | is.na(rep2.matrix[j,k]) | rep1.matrix[j,k]<2 | rep2.matrix[j,k]<2){
#       j+1
#     } else if (rep1.matrix[j,k] > (log2(1.5) + rep2.matrix[j,k]) | (rep2.matrix[j,k] > (log2(1.5) + rep1.matrix[j,k])) == TRUE) 
#     {
#       trans.norm.avg[j,k] <- NA
#     }
#   }
# }
# remove(j,k)
# #check how many NA values there are with above method --> 3.16%
# sum(is.na(trans.norm.avg))/length(c(trans.norm.avg))

#Changing to 2x instead of 1.5x as cutoff for difference between reps
for(k in 1:ncol(rep1.matrix)){
  for(j in 1:nrow(rep1.matrix)){
    if(is.na(rep1.matrix[j,k]) | is.na(rep2.matrix[j,k]) | rep1.matrix[j,k]<2 | rep2.matrix[j,k]<2){
      j+1
    } else if (rep1.matrix[j,k] > (log2(2) + rep2.matrix[j,k]) | (rep2.matrix[j,k] > (log2(2) + rep1.matrix[j,k])) == TRUE) 
    {
      trans.norm.avg[j,k] <- NA
    }
  }
}
remove(j,k)
#check how many NA values there are with above method --> 2.13%
sum(is.na(trans.norm.avg))/length(c(trans.norm.avg))
#what percent of values above the threshold (2) are now NA? --> 3.49%
sum(is.na(c(trans.norm.avg[above2])))/length(above2)
#what about percent if still using 1.5 as the threshold? --> 9.81%
#what about original formula(cutoff 0.2 and 1.5) --> 10.27% NA above 0.2

# For now, going forward with Patrick's rule modified to threshold of 2 for application
# and 2x as the rule (1 on log2)
  
#subtract GST - subtract GST at the same dilution as the other antigen was diluted
#need to deal with NAs this time
  # subtract GST from all except AMA1; that means MSP1-19, Hyp2, GEXP18,

#merge target metadata with trans.norm.avg (transposed)
target.df <- merge(target_meta.df, t(trans.norm.avg), by.x = "Name", by.y = "row.names", sort = FALSE)

#Extract GST elements for all print buffers and concentrations
GST <- c(grep("GST", target.df$Name))
GST.df <- target.df[GST,]

#Plot all GST data
GSTval.df <- GST.df[,(ncol(GST.df) - ncol(t(trans.norm.avg))+1):ncol(GST.df)]
row.names(GSTval.df) <- GST.df$Name

GST_val <- c(as.matrix(GSTval.df))

png(filename = paste0(study, "_GST_ALL.tif"), width = 4, height = 4, units = "in", res = 600)
par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))
plot(GST_val, pch='*', col = "blue", ylim=c(0,max(GST_val,  na.rm = TRUE)*1.25), main = "All GST",
     ylab="Normalized log2(MFI)", xlab="GST averaged duplicates")

#print text on the plots for number of samples where GST and CD4 are above buffer
mtext(paste("Total Samples with GST > 0:", round(sum(GSTval.df > 0, na.rm = TRUE), digits=2), 
            "(", round(sum(GSTval.df > 0, na.rm = TRUE)/length(GST_val)*100, digits=2), "%)"), side=1, cex=0.8, line=0.5, outer=TRUE, xpd=NA, adj=0)

graphics.off()

#function for subtracting GST from a concentration and print buffer
subtractGST <- function(concentration, print_buffer){
  #isolate data for that concentration and print buffer
  con.df <- filter(target.df, Concentration == concentration, Print_Buffer == print_buffer, Expression_Tag == "GST")
  con.df <- tibble::column_to_rownames(con.df, var="Name")
  con.df <- con.df[,(ncol(con.df) - ncol(t(trans.norm.avg))+1):ncol(con.df)]
  
  GSTsub.df <- filter(GST.df, Concentration == concentration, Print_Buffer == print_buffer)
  GSTsub.df <- tibble::column_to_rownames(GSTsub.df, var="Name")
  GSTsub.df <- GSTsub.df[,(ncol(GSTsub.df) - ncol(t(trans.norm.avg))+1):ncol(GSTsub.df)]
  
  con2.df <- data.frame(matrix(0, nrow = nrow(con.df), ncol = ncol(con.df)))
  rownames(con2.df) <- rownames(con.df)
  colnames(con2.df) <- colnames(con.df)
  
  #subtract GST!
  for(b in 1:ncol(con2.df)){
    for(a in 1:nrow(con2.df)){
      #if GST value is NA or main value is NA, then cannot do subtraction and cannot compare with other data
      #so set those values to NA
      if(is.na(GSTsub.df[1,b])|is.na(con.df[a,b])){
        con2.df[a,b] <- NA
        # when the GST value is positive only, subtract GST value, 
        #otherwise want to leave as whatever the value was before (because GST was at or below buffer mean for that sample)
      } else if (GSTsub.df[1,b] > 0){
        #calculate difference in original MFI form (not log2)
        con2.df[a,b] <- 2^con.df[a,b] - 2^GSTsub.df[1,b]
        #can only do log2 if the difference is greater than 0, otherwise set to 0 (normalized log2 value is not above buffer)
        if (con2.df[a,b] > 0) {
          con2.df[a,b] <- log2(con2.df[a,b])
          #if the log2 tag-subtracted value is negative, means normalized value is below buffer mean,
          #so also need to set those negatives to 0 again.
          if(con2.df[a,b] < 0){
            con2.df[a,b] <- 0
          }
        } else { 
          con2.df[a,b] <- 0
        }
      } else {
        con2.df[a,b] <- con.df[a,b]
      }
    }
  }
  remove(a,b)
  
  #return GST-subtracted data
  return(con2.df)
  
}

#Subtract GST from each print buffer and concentration
#I realize it would be easier to use a loop or something here but I was tired
GST100B1.df <- subtractGST(100,1)
GST100B2.df <- subtractGST(100,2)
GST100B3.df <- subtractGST(100,3)

GST25B1.df <- subtractGST(25,1)
GST25B2.df <- subtractGST(25,2)
GST25B3.df <- subtractGST(25,3)

GST6.25B1.df <- subtractGST(6.25,1)
GST6.25B2.df <- subtractGST(6.25,2)
GST6.25B3.df <- subtractGST(6.25,3)

GST1.56B1.df <- subtractGST(1.56,1)
GST1.56B2.df <- subtractGST(1.56,2)
GST1.56B3.df <- subtractGST(1.56,3)

GST0.39B1.df <- subtractGST(0.39,1)
GST0.39B2.df <- subtractGST(0.39,2)
GST0.39B3.df <- subtractGST(0.39,3)

GST0.10B1.df <- subtractGST(0.10, 1)
GST0.10B2.df <- subtractGST(0.10, 2)
GST0.10B3.df <- subtractGST(0.10, 3)

#combine all data back together again
#filter out the GST tagged antigens
no_tags.df <- filter(target.df, !(Expression_Tag == "GST"))
no_tags.df <- tibble::column_to_rownames(no_tags.df, var="Name")
no_tags.df <- no_tags.df[,(ncol(no_tags.df) - ncol(t(trans.norm.avg))+1):ncol(no_tags.df)]
#then rbind the GST data frames to that one. The order of the targets shouldn't matter anymore
optimization.df <- rbind(no_tags.df, GST0.10B1.df, GST0.10B2.df, GST0.10B3.df,
    GST0.39B1.df, GST0.39B2.df, GST0.39B3.df, GST1.56B1.df, GST1.56B2.df, GST1.56B3.df, 
    GST100B1.df, GST100B2.df, GST100B3.df, GST25B1.df, GST25B2.df, GST25B3.df,
    GST6.25B1.df, GST6.25B2.df, GST6.25B3.df)
#Check that dimensions are the same as before
dim(optimization.df) == dim(t(trans.norm.avg))
# View(optimization.df)

# ratio of positive to negative - subtract the negative log2 value for each condition
optimization.df <- t(optimization.df)

Neg <- optimization.df[c(grep("Neg", rownames(optimization.df))),]
CP3 <- optimization.df[c(grep("CP3", rownames(optimization.df))),]
PRISM <- optimization.df[c(grep("PRISM", rownames(optimization.df))),]
Swazi <- optimization.df[c(grep("Swazi", rownames(optimization.df))),]

subtractNeg <- function(x){
for(i in 1:nrow(Neg)){
  for(k in 1:ncol(Neg)){
    if(is.na(Neg[i,k])|is.na(x[i,k])){
      x[i,k] <- NA
    } else {
      x[i,k] <- x[i,k] - Neg[i,k]
    }
  }
}
return(x)
}


#Plot GST values for each sample
CP3.GST <- CP3[,c(grep("GST", colnames(CP3)))]
PRISM.GST <- PRISM[,c(grep("GST", colnames(PRISM)))]
Swazi.GST <- Swazi[,c(grep("GST", colnames(Swazi)))]
Neg.GST <- Neg[,c(grep("GST", colnames(Swazi)))]

png(filename = paste0(study, "_GST_Sample.tif"), width = 7.5, height = 9, units = "in", res = 600)
par(mfrow=c(2,2), oma=c(4,1,1,1),mar=c(6.1,4.1,1.1,1.1))
boxplot(CP3.GST, pch='*', col = "light blue", ylim=c(0,7), main = "CP3 GST",
     ylab="Normalized log2(MFI)", las=2, cex.axis = 0.5)

boxplot(PRISM.GST, pch='*', col = "light blue", ylim=c(0,7), main = "PRISM GST",
     ylab="Normalized log2(MFI)", las=2, cex.axis = 0.5)

boxplot(Swazi.GST, pch='*', col = "light blue", ylim=c(0,7), main = "Swazi GST",
     ylab="Normalized log2(MFI)",las=2, cex.axis = 0.5)

boxplot(Neg.GST, pch='*', col = "light blue", ylim=c(0,7), main = "Neg GST",
     ylab="Normalized log2(MFI)", las=2, cex.axis = 0.5)

#print text on the plots for number of samples where GST and CD4 are above buffer
mtext(paste("GST > 0: CP3 =", round(sum(CP3.GST > 0, na.rm = TRUE), digits=2), 
    "(", round(sum(CP3.GST > 0, na.rm = TRUE)/length(CP3.GST)*100, digits=2), "%), PRISM = ",
    round(sum(PRISM.GST > 0, na.rm = TRUE), digits=2), 
    "(", round(sum(PRISM.GST > 0, na.rm = TRUE)/length(PRISM.GST)*100, digits=2), "%), Swazi =",
    round(sum(Swazi.GST > 0, na.rm = TRUE), digits=2), 
    "(", round(sum(Swazi.GST > 0, na.rm = TRUE)/length(Swazi.GST)*100, digits=2), "%), Neg = ",
    round(sum(Neg.GST > 0, na.rm = TRUE), digits=2), 
    "(", round(sum(Neg.GST > 0, na.rm = TRUE)/length(Neg.GST)*100, digits=2), "%)"), 
    side=1, cex=0.8, line=1.5, outer=TRUE, xpd=NA, adj=0)

graphics.off()

#Now a negative value means the positive is less than the negative control
#Do not set these values to 0
CP3neg <- subtractNeg(CP3)
PRISMneg <- subtractNeg(PRISM)
Swazineg <- subtractNeg(Swazi)

Finaldata <- rbind(CP3neg,PRISMneg,Swazineg)

#get the data in the right format for ggplot2 and linear mixed models

#how to get print buffer as a factor for each antigen
#merge transposed final data with target metadata again
#separate data by print buffer in three data frames
#replace row names with new rownames that are the same for all print buffers, don't include spaces
#bind these back together.
finaltarget.df <- merge(target_meta.df, t(Finaldata), by.y="row.names", by.x = "Name", sort = FALSE)
PB1 <- filter(finaltarget.df, Print_Buffer == "1")
PB2 <- filter(finaltarget.df, Print_Buffer == "2")
PB3 <- filter(finaltarget.df, Print_Buffer == "3")

rownames(PB1) <- PB1$Name
rownames(PB2) <- PB1$Name
rownames(PB3) <- PB1$Name

#Print buffers: 1 = AJ Buffer C		2 = AJ Glycerol buffer		3 = Nexterion Spot
PB1 <- as.data.frame(t(PB1[,6:ncol(PB1)]))
  PB1$print_buffer <- "AJ"
PB2 <- as.data.frame(t(PB2[,6:ncol(PB2)]))
  PB2$print_buffer <- "AJ_Glycerol"
PB3 <- as.data.frame(t(PB3[,6:ncol(PB3)]))
  PB3$print_buffer <- "Nexterion_Spot"

#merge these with sample_meta_f to get other columns, then rbind all together
PB1_meta <- merge(sample_meta_f.df, PB1, by.x = "sample_id", by.y = "row.names",sort=FALSE)
PB2_meta <- merge(sample_meta_f.df, PB2, by.x = "sample_id", by.y = "row.names",sort=FALSE)
PB3_meta <- merge(sample_meta_f.df, PB3, by.x = "sample_id", by.y = "row.names",sort=FALSE)

AHHHH.df <- rbind(PB1_meta, PB2_meta, PB3_meta)

#need to change the sample_id to only CP3, PRISM, and Swazi 
#so that the model knows they are the same sample
AHHHH.df$sample <- "CP3"
AHHHH.df$sample[c(grep("PRISM", AHHHH.df$sample_id))] <- "PRISM"
AHHHH.df$sample[c(grep("Swazi", AHHHH.df$sample_id))] <- "Swazi"

#the names of the columns have spaces. So they cannot be used in lmer and other functions
newnames <- make.names(colnames(AHHHH.df))
colnames(AHHHH.df) <- newnames

#####Plotting#######

#Plot the data for the PRISM and CP3 samples only, organized by geometric mean

AHHsub <- filter(AHHHH.df, sample == "CP3" | sample == "PRISM")
CP3all <- filter(AHHHH.df, sample == "CP3")
CP3sub <- CP3all[11:46]
PRISMall <- filter(AHHHH.df, sample == "PRISM")
PRISMsub <- PRISMall[11:46]

geomean <- (CP3sub + PRISMsub)/2
geomean <- tibble:: rowid_to_column(geomean)
PRISMall <- tibble::rowid_to_column(PRISMall)
CP3all <- tibble::rowid_to_column(CP3all)
AHHsub2 <- rbind(PRISMall, CP3all)

#plot data from top 48 geomeans (10%) of conditions
geomeanAMA1.100 <- geomean[order(geomean$X13_1.PfAMA1.100ug.ml_1, decreasing = TRUE),]
rowidorder <- c(as.character(geomeanAMA1.100$rowid))
geomeanAMA1.100$rowid <- factor(geomeanAMA1.100$rowid, levels = unique(rowidorder))

AHHsub2$rowid <- factor(AHHsub2$rowid, levels = unique(rowidorder))

#export table of conditions matching top 10% rowid. 

png(filename = paste0("AMA1.100.GM.tif"), width = 8, height = 3, units = "in", res = 1200)
par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))

ggplot(geomeanAMA1.100[1:48,], aes(x = rowid, y = X13_1.PfAMA1.100ug.ml_1)) + theme_bw() + 
  geom_point() + ylab("Normalized Log2(MFI)") + theme(axis.text.x = element_text( vjust = 1, size = 10))

graphics.off()

png(filename = paste0("AMA1.100.2samples.tif"), width = 8, height = 3, units = "in", res = 1200)
par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))

ggplot(AHHsub2, aes(x = rowid, y = X13_1.PfAMA1.100ug.ml_1, color = sample)) + theme_bw() + 
  geom_point() + ylab("Normalized Log2(MFI)") + theme(axis.text.x = element_text( vjust = 1, size = 10))

graphics.off()






#violin, scatter, or boxplots of the data

#all samples and all variables - AMA1 100 - Scatter plot
png(filename = paste0(study, "AMA1.100_facet_grid.tif"), width = 8, height = 5, units = "in", res = 1200)
par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))

ggplot(AHHHH.df, aes(x=slide_type, y=X13_1.PfAMA1.100ug.ml_1, color = sample, shape = blocking_buffer)) + 
  geom_point(size = 1)  + 
  theme_bw() + 
  theme(text = element_text(size=9), axis.text.x = element_text(angle = 90, hjust=0.95, vjust = 0.1)) + 
  labs(x = "Slide Type", y = "Normalized Log2(Positive/Negative)", title = "AMA1 (100 µg/mL)") +
  facet_grid(print_buffer ~ block_dilution )

graphics.off()

#violin plot of slide type and print buffer, everything else combined 
png(filename = paste0(study, "AMA1.100_ST_PB.tif"), width = 8, height = 5, units = "in", res = 1200)
par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))

ggplot(AHHHH.df, aes(x=reorder(slide_type, X13_1.PfAMA1.100ug.ml_1, FUN=max), y=X13_1.PfAMA1.100ug.ml_1, color = print_buffer, fill = print_buffer)) + 
  geom_violin(trim = FALSE) + 
  theme_bw() +
  guides(color=FALSE) +
  theme(text = element_text(size=9), axis.text.x = element_text(angle = 45, hjust = 1), legend.position="top") +
  labs(x = "Slide Type", y = "Normalized Log2(Positive/Negative)", title = "AMA1 (100 µg/mL)", fill = "Print Buffer")

graphics.off()

#Make a function of the plots to use with other columns

#column is the column name from AHHHH.df for each antigen and dilution
#antigen and dilution are character vectors 
# antigenplots <- function(column, antigen, dilution){
#   
# #all samples and all variables - Scatter plot faceted
# png(filename = paste0(antigen,"_", dilution, "_facet_grid.tif"), width = 8, height = 5, units = "in", res = 1200)
# par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))
# 
# p1 <- ggplot(AHHHH.df, aes(x=slide_type, y=column, color = sample, shape = blocking_buffer)) + 
#   geom_point(size = 1)  + 
#   theme_bw() + 
#   theme(text = element_text(size=9), axis.text.x = element_text(angle = 90, hjust=0.95, vjust = 0.1)) + 
#   labs(x = "Slide Type", y = "Normalized Log2(Positive/Negative)", title = paste(antigen, dilution, "µg/mL")) +
#   facet_grid(print_buffer ~ block_dilution )
# 
# print(p1)
# graphics.off()
# 
# #violin plot of slide type and print buffer, everything else combined 
# png(filename = paste0(antigen,"_", dilution, "_ST_PB.tif"), width = 8, height = 5, units = "in", res = 1200)
# par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))
# 
# p2 <- ggplot(AHHHH.df, aes(x=slide_type, y=column, color = print_buffer, fill = print_buffer)) + 
#   geom_violin(trim = FALSE) + 
#   theme_bw() +
#   guides(color=FALSE) +
#   theme(text = element_text(size=9), axis.text.x = element_text(angle = 45, hjust = 1), legend.position="top") +
#   labs(x = "Slide Type", y = "Normalized Log2(Positive/Negative)", title = paste(antigen, dilution, "µg/mL"), fill = "Print Buffer")
# 
# print(p2)
# graphics.off()
# 
# }
# 
# antigenplots(X13_1.PfAMA1.100ug.ml_1, "AMA1", "100")

#heatmap of all the data 
#want a heatmap with antigen name and dilution on the top (x)
#and the combo on the side. Sorted by highest total intensity (sum of the rows)

#prepare relevant data as a matrix
heatdata <- as.matrix(AHHHH.df[,11:(ncol(AHHHH.df)-2)])

#prepare vector for order of rows and columns
#AHHHH.df is currently ordered by sample > print buffer > slide type > blocking buffer > dilution
#and the columns are ordered by antigen > dilution

#prepare character vector labels? or make labels later?
png(filename = paste0(study, "heatmap2.All.tif"), width = 4.5, height = 8, units = "in", res = 1200)
par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))

heatmap1(heatdata, row.order = NULL)

graphics.off()

#heatmap of the data separated by print buffer 
AJ <- filter(AHHHH.df, print_buffer == "AJ")
AJ.mat <- as.matrix(AJ[,11:46])

AJGly <- filter(AHHHH.df, print_buffer == "AJ_Glycerol")
AJGly.mat <- as.matrix(AJGly[,11:46])

Next <- filter(AHHHH.df, print_buffer == "Nexterion_Spot")
Next.mat <- as.matrix(Next[,11:46])

png(filename = paste0(study, "AJ_heatmap.All.tif"), width = 4.5, height = 8, units = "in", res = 1200)
par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))
heatmap1(AJ.mat, row.order = NULL)
graphics.off()



  
# 
# heatmap(as.matrix(data))
# 
# #heatmap with ggplot2 - not done yet
# ggplot(AHHHH.df, aes(slide_type, blocking_buffer, block_dilution)) +
#   geom_tile(aes(fill = AHHHH.df[11]), color = "white") +
#   scale_fill_gradient(low = "red", high = "steelblue") +
#   ylab("") +
#   xlab("") +
#   theme(legend.title = element_text(size = 9),
#         legend.text = element_text(size = 9),
#         plot.title = element_text(size=9),
#         axis.title=element_text(size=9,face="bold"),
#         axis.text.x = element_text(angle = 90, hjust = 1)) +
#   labs(fill = "Normalized Log2(MFI)")
# 
# #hierarchical clustering - this isn't working
# clusters <- hclust(dist(AHHHH.df[11]), na.rm = TRUE)
# plot(clusters)

######### Linear Mixed Effects Models ########

#once the data is formatted right, do linear mixed model for each antigen
#fixed effects are: slide_type, print_buffer, blocking_buffer, block_dilution
#random effects are: sample; using random slope model where print buffer and blocking buffer (and dilution) affect sample variation

#save workspace image to start from this point
save.image(file = "OptimizationLMMready.RData")

#AMA1, 100 ug/mL, - change REML to FALSE for anova
colnames(AHHHH.df[11])
fullmodel <- lmer(X13_1.PfAMA1.100ug.ml_1 ~ slide_type + blocking_buffer 
    + block_dilution + print_buffer + (1|sample), 
    REML = FALSE, data = AHHHH.df)
summary(fullmodel)
r.squaredGLMM(fullmodel)

#look at the impact of each individual factor by removing 1 factor at a time
#print buffer
model_PB <- lmer(X13_1.PfAMA1.100ug.ml_1 ~ slide_type + blocking_buffer 
               + block_dilution + (1|sample), 
               REML = FALSE, data = AHHHH.df)
  summary(model_PB)
  r.squaredGLMM(model_PB)
  #compare these with likelihood ratio test
  anova(model_PB, fullmodel)
  
#slide type 
  model_ST <- lmer(X13_1.PfAMA1.100ug.ml_1 ~ print_buffer + blocking_buffer 
                   + block_dilution + (1|sample), 
                   REML = FALSE, data = AHHHH.df)
  summary(model_ST)
  r.squaredGLMM(model_ST)
  #compare these with likelihood ratio test
  anova(model_ST, fullmodel)  
  
#blocking buffer
  model_BB <- lmer(X13_1.PfAMA1.100ug.ml_1 ~ print_buffer + slide_type 
                   + block_dilution + (1|sample), 
                   REML = FALSE, data = AHHHH.df)
  summary(model_BB)
  r.squaredGLMM(model_BB)
  #compare these with likelihood ratio test
  anova(model_BB, fullmodel)
  
#block dilution
  model_BD <- lmer(X13_1.PfAMA1.100ug.ml_1 ~ print_buffer + slide_type 
                  + blocking_buffer + (1|sample), 
                   REML = FALSE, data = AHHHH.df)
  summary(model_BD)
  r.squaredGLMM(model_BD)
  #compare these with likelihood ratio test
  anova(model_BD, fullmodel)
  
#look at the 2 way interactions between factors 
  int2model <- lmer(X13_1.PfAMA1.100ug.ml_1 ~ slide_type + blocking_buffer + block_dilution + print_buffer +
                      (slide_type + blocking_buffer + block_dilution + print_buffer)^2 +
                      (1|sample), REML = FALSE, data = AHHHH.df)
  summary(int2model)
  r.squaredGLMM(int2model)
  
  anova(int2model)
  
  #compare these with likelihood ratio test
  anova(fullmodel, int2model)
  
  #Plot residuals - int2model
  png(filename = "AMA1.100.Int2.Res.tif", width = 8, height = 4.5, units = "in", res = 1200)
  par(mfrow=c(1,2), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))
  
  plot(fitted(int2model),residuals(int2model),  pch='*', col = "blue", xlab = "Fitted", ylab = "Residuals")
  abline(a=0, b=0)
  hist(residuals(int2model), xlab = "Residuals")
  
  graphics.off()
  
#look at 3 way interactions between factors
  int3model <- lmer(X13_1.PfAMA1.100ug.ml_1 ~ slide_type + blocking_buffer + block_dilution + print_buffer +
                      (slide_type + blocking_buffer + block_dilution + print_buffer)^3 +
                      (1|sample), REML = FALSE, data = AHHHH.df)
  summary(int3model)
  r.squaredGLMM(int3model)
  
  anova(int3model)
  
  #compare these with likelihood ratio test
  anova(fullmodel, int3model)
  anova(int2model, int3model)
  
  #Plot residuals - int3model
  png(filename = "AMA1.100.Int3.Res.tif", width = 8, height = 4.5, units = "in", res = 1200)
  par(mfrow=c(1,2), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))
  
  plot(fitted(int3model),residuals(int3model),  pch='*', col = "blue", xlab = "Fitted", ylab = "Residuals")
  abline(a=0, b=0)
  hist(residuals(int3model), xlab = "Residuals")
  
  graphics.off()
  
#look at 4 way interactions between factors
  int4model <- lmer(X13_1.PfAMA1.100ug.ml_1 ~ slide_type + blocking_buffer + block_dilution + print_buffer +
                      (slide_type + blocking_buffer + block_dilution + print_buffer)^4 +
                      (1|sample), REML = FALSE, data = AHHHH.df)
  summary(int4model)
  r.squaredGLMM(int4model)
  
  anova(int4model)
  
  #compare these with likelihood ratio test
  anova(fullmodel, int4model)
  anova(int2model, int4model)
  anova(int3model, int4model)
  
  #Plot residuals - int4model
  png(filename = "AMA1.100.Int4.Res.tif", width = 8, height = 4.5, units = "in", res = 1200)
  par(mfrow=c(1,2), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))
  
  plot(fitted(int4model),residuals(int4model),  pch='*', col = "blue", xlab = "Fitted", ylab = "Residuals")
  abline(a=0, b=0)
  hist(residuals(int4model), xlab = "Residuals")
  
  graphics.off()
  
#modified 3 way interaction model - significant 3 way interactions taken from int3model
  int3.2 <- lmer(X13_1.PfAMA1.100ug.ml_1 ~ slide_type + blocking_buffer + block_dilution + print_buffer +
                   slide_type*block_dilution*print_buffer + slide_type*blocking_buffer*block_dilution  + 
                   (1|sample), REML = FALSE, data = AHHHH.df)

  summary(int3.2)
  r.squaredGLMM(int3.2)
  anova(int3.2)
  
  #compare these with likelihood ratio test
  anova(int2model, int3.2)
  anova(int3model, int3.2)
  
  #Plot residuals - int3.2
  png(filename = "AMA1.100.Int3.2.Res.tif", width = 8, height = 4.5, units = "in", res = 1200)
  par(mfrow=c(1,2), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))
  
  plot(fitted(int3.2),residuals(int3.2),  pch='*', col = "blue", xlab = "Fitted", ylab = "Residuals")
  abline(a=0, b=0)
  hist(residuals(int3.2), xlab = "Residuals")
  
  graphics.off()
  
  
#prepare a better organized summary of the model and export to a file.
#summary can only be written to a text file, and doesn't keep columns organized
tidy(fullmodel)
augment(fullmodel)

write.csv(tidy(fullmodel), file = "AMA1.100.LMERTidy.csv")
write.csv(augment(fullmodel), file = "AMA1.100.LMERAug.csv")

#test for main effects and interactions with Likelihood Ratio Test?

intmodel <- lmer(X13_1.PfAMA1.100ug.ml_1 ~ slide_type * blocking_buffer 
                    * block_dilution * print_buffer + (1|sample), 
                    REML = FALSE, data = AHHHH.df)
summary(intmodel)

#Plot residuals
png(filename = "AMA1.100.Residuals.tif", width = 8, height = 4.5, units = "in", res = 1200)
par(mfrow=c(1,2), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))

plot(fitted(intmodel),residuals(intmodel),  pch='*', col = "blue", xlab = "Fitted", ylab = "Residuals")
abline(a=0, b=0)
hist(residuals(intmodel), xlab = "Residuals")

graphics.off()


r.squaredGLMM(intmodel)

write.csv(tidy(intmodel), file = "AMA1.100.LMERTidyINT.csv")
write.csv(augment(intmodel), file = "AMA1.100.LMERAugINT.csv")

anova(fullmodel, intmodel)

anova(intmodel)

#get significance and do pairwise comparisons using emmeans
#default adjustment for multiple comparisons is Tukey
#main effects first, however we are not interested in main effects :(
emmeans(fullmodel, pairwise ~ print_buffer)
emmeans(fullmodel, pairwise ~ slide_type)
emmeans(fullmodel, pairwise ~ blocking_buffer)
emmeans(fullmodel, pairwise ~ block_dilution)

#calculate emms
emm.intmodel <- emmeans(intmodel, ~ print_buffer | slide_type | blocking_buffer | block_dilution)

#repeat for int3.2
emm3.2 <- emmeans(int3.2, ~ print_buffer | slide_type | blocking_buffer | block_dilution)

#plot the data
png(filename = "AMA1.100.emmip.tif", width = 8, height = 5, units = "in", res = 1200)
par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))
emmip(intmodel, print_buffer ~ block_dilution | blocking_buffer | slide_type) + element_text(size=8)

graphics.off()

#plot emmeans for int3.2 model
png(filename = "AMA1.100.emmip3.2.tif", width = 8, height = 5, units = "in", res = 1200)
par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))
emmip(emm3.2, print_buffer ~ block_dilution | blocking_buffer | slide_type)

graphics.off()

#extract the emmeans and sort highest to lowest
emm3.2.df <- as.data.frame(emm3.2)
emm3.2.df <- emm3.2.df[order(emm3.2.df$emmean, decreasing = TRUE),]

#make rowid a column to match with AHHHH.df
emm3.2.df <- tibble::rowid_to_column(emm3.2.df)

#merge with AHHHH.df
emmAMA1.100.df <- merge(emm3.2.df, AHHHH.df, all.x = TRUE, sort = FALSE)

#keep same order that we have in the excel file (sorted highest to lowest EMM)
#maybe try to use rownumber as factor, sorted by EMM? and then put in better labels later?
rowidorder <- c(as.character(emmAMA1.100.df$rowid))
emmAMA1.100.df$rowid <- factor(emmAMA1.100.df$rowid, levels = rev(unique(rowidorder)))

#plot original data from groups with higest emmeans (top 36)
emmAMA1.100.df$rowid <- factor(emmAMA1.100.df$rowid, levels = unique(rowidorder))

emmAMA1.100sub <- emmAMA1.100.df[1:108,]

png(filename = paste0("H.AMA1.100.3.2.tif"), width = 8, height = 3, units = "in", res = 1200)
par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))

ggplot(emmAMA1.100sub, aes(x = rowid, y = X13_1.PfAMA1.100ug.ml_1, color = sample)) + theme_bw() + geom_point() + ylab("Normalized Log2(MFI)") + theme(axis.text.x = element_text( vjust = 1, size = 10))

graphics.off()




#get letters for all pairwise comparisons, interaction model

#this took several hours (>4) to run, and eventually finished. there are too many 
#combinations (480 pairwise), but it worked!!! 
#the contrast function used by cld function automatically uses Tukey adjustment
#for multiple comparisons
pairwiseletters <- cld(emm.intmodel, by = NULL, Letters = LETTERS, sort = TRUE, reversed = TRUE, details = TRUE)
#save the results to a file. 
write.csv(as.data.frame(pairwiseletters$emmeans), file = "OptimizationALLPairwiseEmmeans.csv")
write.csv(as.data.frame(pairwiseletters$comparisons), file = "OptimizationALLPairwiseComparisons.csv")
#save the workspace
save.image("OptimizationAfterPairwise.RData")

#MSP1-19, 100 ug/mL, 
colnames(AHHHH.df[17])
fullmodel <- lmer(X19_1.PfMSP1.19.100ug.ml_1 ~ slide_type + blocking_buffer 
                  + block_dilution + print_buffer + (1|sample), 
                  REML = FALSE, data = AHHHH.df)
summary(fullmodel)
r.squaredGLMM(fullmodel)

#plot residuals - they look mostly good, some heteroskedasticity though and data looks linear
plot(fitted(fullmodel),residuals(fullmodel))
#check normality - looks good
hist(residuals(fullmodel))

#prepare a better organized summary of the model and export to a file.
#summary can only be written to a text file, and doesn't keep columns organized
write.csv(tidy(fullmodel), file = "MSP1-19.100.LMERTidy.csv")
write.csv(augment(fullmodel), file = "MSP1-19.100.LMERAug.csv")

#test for main effects and interactions with Likelihood Ratio Test?
intmodel <- lmer(X19_1.PfMSP1.19.100ug.ml_1 ~ slide_type * blocking_buffer 
                 * block_dilution * print_buffer + (1|sample), 
                 REML = FALSE, data = AHHHH.df)
summary(intmodel)

#Plot residuals 
png(filename = "MSP1.19.100.Residuals.tif", width = 7.8, height = 4.5, units = "in", res = 1200)
par(mfrow=c(1,2), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))

plot(fitted(intmodel),residuals(intmodel),  pch='*', col = "blue", xlab = "Fitted", ylab = "Residuals")
abline(a=0, b=0)
hist(residuals(intmodel), xlab = "Residuals")

graphics.off()

r.squaredGLMM(intmodel)

write.csv(tidy(intmodel), file = "MSP1.19.100.LMERTidyINT.csv")
write.csv(augment(intmodel), file = "MSP1.19.100.LMERAugINT.csv")

anova(fullmodel, intmodel)

anova(intmodel)

#get significance and do pairwise comparisons using emmeans
#default adjustment for multiple comparisons is Tukey
#main effects first, however we are not interested in main effects :(
emmeans(fullmodel, pairwise ~ print_buffer)
emmeans(fullmodel, pairwise ~ slide_type)
emmeans(fullmodel, pairwise ~ blocking_buffer)
emmeans(fullmodel, pairwise ~ block_dilution)

#plot the data
png(filename = "MSP1-19.100.emmip.tif", width = 8, height = 5, units = "in", res = 1200)
par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))
emmip(intmodel, print_buffer ~ block_dilution | blocking_buffer | slide_type)

graphics.off()

#get letters for all pairwise comparisons, interaction model
emm.intmodel <- emmeans(intmodel, ~ print_buffer | slide_type | blocking_buffer | block_dilution)

#this took several hours (>4) to run, and eventually finished. there are too many 
#combinations (480 pairwise), but it worked!!! 
#the contrast function used by cld function automatically uses Tukey adjustment
#for multiple comparisons
pairwiseletters <- cld(emm.intmodel, by = NULL, Letters = LETTERS, sort = TRUE, reversed = TRUE, details = TRUE)
#save the results to a file. 
write.csv(as.data.frame(pairwiseletters$emmeans), file = "MSP1-19.100.ALLPairwiseEmmeans.csv")
write.csv(as.data.frame(pairwiseletters$comparisons), file = "MSP1-19.100.ALLPairwiseComparisons.csv")
#above, forgot to change the name of pairwiseletters to be antigen specific!
#so this pairwise overwrote the previous one for AMA1

#save the workspace
save.image("OptimizationAfterPairwiseMSP1-19.RData")

#Hyp2, 100 ug/mL, 
colnames(AHHHH.df[23])
fullmodel <- lmer(X25_1.Hyp2.100ug.ml_1 ~ slide_type + blocking_buffer 
                  + block_dilution + print_buffer + (1|sample), 
                  REML = FALSE, data = AHHHH.df)
summary(fullmodel)
r.squaredGLMM(fullmodel)

#plot residuals - they look mostly good, some heteroskedasticity though and data looks linear
plot(fitted(fullmodel),residuals(fullmodel))
#check normality - looks good
hist(residuals(fullmodel))

#prepare a better organized summary of the model and export to a file.
#summary can only be written to a text file, and doesn't keep columns organized
write.csv(tidy(fullmodel), file = "Hyp2.100.LMERTidy.csv")
write.csv(augment(fullmodel), file = "Hyp2.100.LMERAug.csv")

#test for main effects and interactions with Likelihood Ratio Test?
intmodel <- lmer(X25_1.Hyp2.100ug.ml_1 ~ slide_type * blocking_buffer 
                 * block_dilution * print_buffer + (1|sample), 
                 REML = FALSE, data = AHHHH.df)
summary(intmodel)

#Plot residuals 
png(filename = "Hyp2.Residuals.tif", width = 7.8, height = 4.5, units = "in", res = 1200)
par(mfrow=c(1,2), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))

plot(fitted(intmodel),residuals(intmodel),  pch='*', col = "blue", xlab = "Fitted", ylab = "Residuals")
abline(a=0, b=0)
hist(residuals(intmodel), xlab = "Residuals")

graphics.off()

r.squaredGLMM(intmodel)

write.csv(tidy(intmodel), file = "Hyp2.100.LMERTidyINT.csv")
write.csv(augment(intmodel), file = "Hyp2.100.LMERAugINT.csv")

anova(fullmodel, intmodel)

anova(intmodel)

#plot the model data
png(filename = "Hyp2.100.emmip.tif", width = 8, height = 5, units = "in", res = 1200)
par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))
emmip(intmodel, print_buffer ~ block_dilution | blocking_buffer | slide_type)

graphics.off()

#get letters for all pairwise comparisons, interaction model
emm.intmodel <- emmeans(intmodel, ~ print_buffer | slide_type | blocking_buffer | block_dilution)

#this took several hours (>6.5) to run, and eventually finished. there are too many 
#combinations (480 pairwise), but it worked!!! 
#the contrast function used by cld function automatically uses Tukey adjustment
#for multiple comparisons
Hyp2pairwiseletters <- cld(emm.intmodel, by = NULL, Letters = LETTERS, sort = TRUE, reversed = TRUE, details = TRUE)
#save the results to a file. 
write.csv(as.data.frame(Hyp2pairwiseletters$emmeans), file = "Hyp2.100.ALLPairwiseEmmeans.csv")
write.csv(as.data.frame(Hyp2pairwiseletters$comparisons), file = "Hyp2.100.ALLPairwiseComparisons.csv")


######### Plotting the data from the successful combos in pairwiseletters #######

#import the data which I decided was good in excel
AMA1.100.best <- read.csv("AMA1.100.BEST.csv")
MSP1.19.100.best <- read.csv("MSP1-19.100.BEST.csv")

#add rowid as a column to get a number for each combination
AMA1.100.best <- tibble::rowid_to_column(AMA1.100.best)
MSP1.19.100.best <- tibble::rowid_to_column(MSP1.19.100.best)

#merge with the AHHHH.df to get original values for all three samples
#some are NA! I didn't realize that samples with NA would be included
AMA1.100.exp <- merge(AMA1.100.best, AHHHH.df, all.x = TRUE, sort = FALSE)
MSP1.19.100.exp <- merge(MSP1.19.100.best, AHHHH.df, all.x = TRUE, sort = FALSE)

#keep same order that we have in the excel file (sorted highest to lowest EMM)
#maybe try to use rownumber as factor, sorted by EMM? and then put in better labels later?
rowidorder <- c(as.character(AMA1.100.exp$rowid))
AMA1.100.exp$rowid <- factor(AMA1.100.exp$rowid, levels = rev(unique(rowidorder)))

rowidorderMSP <- c(as.character(MSP1.19.100.exp$rowid))
MSP1.19.100.exp$rowid <- factor(MSP1.19.100.exp$rowid, levels = rev(unique(rowidorderMSP)))

#COORD_FLIP!!! Plot every point for each combo, different colors for each sample, 
#AMA1 100
png(filename = paste0("AMA1.100.best.dot.tif"), width = 3, height = 8, units = "in", res = 1200)
par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))

ggplot(AMA1.100.exp, aes(x = rowid, y = X13_1.PfAMA1.100ug.ml_1, color = sample)) + theme_bw() + geom_point() + coord_flip() + ylab("Normalized Log2(MFI)") + theme(axis.text.x = element_text( vjust = 1, hjust=1, size = 10))

graphics.off()

#MSP1.19 100
png(filename = paste0("MSP1_19.100.best.dot.tif"), width = 3, height = 8, units = "in", res = 1200)
par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))

ggplot(MSP1.19.100.exp, aes(x = rowid, y = X19_1.PfMSP1.19.100ug.ml_1, color = sample)) + theme_bw() + geom_point() + coord_flip() + ylab("Normalized Log2(MFI)") + theme(axis.text.x = element_text( vjust = 1, hjust=1, size = 10))

graphics.off()

#Same plots - NOT COORD_FLIP to put on PPT slides 
#AMA1 100
AMA1.100.exp$rowid <- factor(AMA1.100.exp$rowid, levels = unique(rowidorder))

png(filename = paste0("H.AMA1.100.best.dot.tif"), width = 8, height = 3, units = "in", res = 1200)
par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))

ggplot(AMA1.100.exp, aes(x = rowid, y = X13_1.PfAMA1.100ug.ml_1, color = sample)) + theme_bw() + geom_point() + ylab("Normalized Log2(MFI)") + theme(axis.text.x = element_text( vjust = 1, size = 10))

graphics.off()

#MSP1.19 100
MSP1.19.100.exp$rowid <- factor(MSP1.19.100.exp$rowid, levels = unique(rowidorderMSP))

png(filename = paste0("H.MSP1_19.100.best.dot.tif"), width = 8, height = 3, units = "in", res = 1200)
par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))

ggplot(MSP1.19.100.exp, aes(x = rowid, y = X19_1.PfMSP1.19.100ug.ml_1, color = sample)) + theme_bw() + geom_point() + ylab("Normalized Log2(MFI)") + theme(axis.text.x = element_text(vjust = 1, size = 10))

graphics.off()

#Which combinations are present in multiple antigens? 
AMA1.100.Combo <- AMA1.100.best[,2:5]
MSP1.19.100.Combo <- MSP1.19.100.best[,2:5]

ComboMult <- merge(AMA1.100.Combo, MSP1.19.100.Combo, sort = FALSE)
#12 combos are present in AMA1.100 and MSP1.19.100 best lists