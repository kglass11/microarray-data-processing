#optimization data analysis 
#KG

setwd("/Users/Katie/Desktop/R files from work/270717 Optimisation/")
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

library(broom)
library(lme4)
library(multcomp)
library(lsmeans)
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

# load("ATT93447.RData")
load("OptimizationAfterPairwise.RData")

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

### Check for deviant technical replicates, automatically exclude (set to NA)
# Use Patrick's formula for ELISA to compare replicates within one array
# if rep1 or rep2 is more than 1.5 times rep2 or rep1, respectively, exclude that pair
# Added: if either rep1 or rep2 is less than 0.2, then don't apply the rule
for(k in 1:ncol(rep1.matrix)){
   for(j in 1:nrow(rep1.matrix)){
     if(is.na(rep1.matrix[j,k]) | is.na(rep2.matrix[j,k]) | rep1.matrix[j,k]<0.2 | rep2.matrix[j,k]<0.2){
       j+1
     } else if (rep1.matrix[j,k] > (log2(1.5) + rep2.matrix[j,k]) | (rep2.matrix[j,k] > (log2(1.5) + rep1.matrix[j,k])) == TRUE) 
     {
       trans.norm.avg[j,k] <- NA
     }
   }
 }
remove(j,k)
  
#subtract GST - subtract GST at the same dilution as the other antigen was diluted
#need to deal with NAs this time
  # subtract GST from all except AMA1; that means MSP1-19, Hyp2, GEXP18,

#merge target metadata with trans.norm.avg (transposed)
target.df <- merge(target_meta.df, t(trans.norm.avg), by.x = "Name", by.y = "row.names", sort = FALSE)

#Extract GST elements for all print buffers and concentrations
GST <- c(grep("GST", target.df$Name))
GST.df <- target.df[GST,]

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
View(optimization.df)

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

#heatmap of all the data - this looks terrible and still has unwanted columns
#want a heatmap of all of the data, with antigen name and dilution on the top (x)
#and the combo on the side. Sorted by highest total intensity (sum of the rows)
data <- AHHHH.df[,sapply(AHHHH.df, is.numeric)]

heatmap(as.matrix(data))

#heatmap with ggplot2 - not done yet
ggplot(AHHHH.df, aes(slide_type, blocking_buffer, block_dilution)) +
  geom_tile(aes(fill = AHHHH.df[11]), color = "white") +
  scale_fill_gradient(low = "red", high = "steelblue") +
  ylab("List of genes ") +
  xlab("List of patients") +
  theme(legend.title = element_text(size = 9),
        legend.text = element_text(size = 9),
        plot.title = element_text(size=9),
        axis.title=element_text(size=9,face="bold"),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(fill = "Normalized Log2(MFI)")

#hierarchical clustering - this isn't working
clusters <- hclust(dist(AHHHH.df[11]), na.rm = TRUE)
plot(clusters)

######### Linear Mixed Effects Models ########

#once the data is formatted right, do linear mixed model for each antigen
#fixed effects are: slide_type, print_buffer, blocking_buffer, block_dilution
#random effects are: sample; using random slope model where print buffer and blocking buffer (and dilution) affect sample variation

#AMA1, 100 ug/mL, 
colnames(AHHHH.df[11])
fullmodel <- lmer(X13_1.PfAMA1.100ug.ml_1 ~ slide_type + blocking_buffer 
    + block_dilution + print_buffer + (1|sample), 
    REML = FALSE, data = AHHHH.df)
summary(fullmodel)
r.squaredGLMM(fullmodel)

#plot residuals - they look good, not heteroskedastic and data looks linear
plot(fitted(fullmodel),residuals(fullmodel))
#check normality - looks good
hist(residuals(fullmodel))

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
plot(fitted(intmodel),residuals(intmodel))
hist(residuals(intmodel))
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

#plot the data
png(filename = "AMA1.100.emmip.tif", width = 8, height = 5, units = "in", res = 1200)
par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))
emmip(intmodel, print_buffer ~ block_dilution | blocking_buffer | slide_type) + element_text(size=8)

graphics.off()

#get letters for all pairwise comparisons, interaction model
emm.intmodel <- emmeans(intmodel, ~ print_buffer | slide_type | blocking_buffer | block_dilution)

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
plot(fitted(intmodel),residuals(intmodel))
hist(residuals(intmodel))
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
