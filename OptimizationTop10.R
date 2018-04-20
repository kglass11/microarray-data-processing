#continuing on from optimizationanalysis.R

#Top 10% conditions based on signal from CP3 and PRISM separately,
  #once with and once without the ratio from positive to negative
  #all will have GST subtraction the same way

conditions <- CP3all[,c(1,9,10,11,48)]

#1. CP3 and PRISM with ratio from positive to negative (CP3all has everything for CP3 samples only)

  #AMA1, 100
  
  #sort/order by desired column, highest to lowest values
  CP3.AMA1.100 <- CP3all[order(CP3all$X13_1.PfAMA1.100ug.ml_1, decreasing = TRUE),]
  rowidorder <- c(as.character(CP3.AMA1.100$rowid))
  CP3.AMA1.100$rowid <- factor(CP3.AMA1.100$rowid, levels = unique(rowidorder))
  
  PRISM.AMA1.100 <- PRISMall[order(PRISMall$X13_1.PfAMA1.100ug.ml_1, decreasing = TRUE),]
  rowidorder <- c(as.character(PRISM.AMA1.100$rowid))
  PRISM.AMA1.100$rowid <- factor(PRISM.AMA1.100$rowid, levels = unique(rowidorder))

  #find top 48 (10%) of conditions and export table of conditions matching top 10% rowid. 
  conditionsAMA1.100 <- merge(CP3.AMA1.100, conditions, sort = FALSE)
  conditionAMA1sub <- conditionsAMA1.100[1:48, 1:5]
  write.csv(conditionAMA1sub, file = "CP3top10AMA1.100.csv")
  
  PRISM.cond.AMA1.100 <- merge(PRISM.AMA1.100, conditions, sort = FALSE)
  PRISM.cond.AMA1sub <- PRISM.cond.AMA1.100[1:48, 1:5]
  write.csv(PRISM.cond.AMA1sub, file = "PRISMtop10AMA1.100.csv")

  #make one data frame with PRISM and CP3 best conditions
  BothAMA1.100 <- rbind(conditionsAMA1.100[1:48,], PRISM.cond.AMA1.100[1:48,])
  
  #plot the top 10% values and rowids for both CP3 and PRISM on one plot - this isn't working yet
  png(filename = paste0("Both.Top.AMA1.100.tif"), width = 8, height = 3, units = "in", res = 1200)
  par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))
  
  base <- ggplot(CP3.AMA1.100[1:48,]) + geom_point(aes(x = rowid, y = X13_1.PfAMA1.100ug.ml_1), col = "red") + theme_bw() + 
    labs(x = "Row ID", y = "Normalized Log2(Positive/Negative)", title = "Best Conditions AMA1 100 µg/mL") + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 9)) +
    ylim(0,9)
  
  base + geom_point(data = PRISM.AMA1.100[1:48,], aes(x = rowid, y = X13_1.PfAMA1.100ug.ml_1), col= "green")
  
  graphics.off()
  
  #plot the top 10% values and rowids for CP3 alone
  png(filename = paste0("CP3.Top.AMA1.100.tif"), width = 8, height = 3, units = "in", res = 1200)
  par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))

  ggplot(CP3.AMA1.100[1:48,], aes(x = rowid, y = X13_1.PfAMA1.100ug.ml_1)) + geom_point(col = "red") + theme_bw() + 
  labs(x = "Row ID", y = "Normalized Log2(Positive/Negative)", title = "CP3 AMA1 100 µg/mL") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 9)) +
  ylim(0,9)

  graphics.off()
  
  #plot the top 10% values and rowids for PRISM alone
  png(filename = paste0("PRISM.Top.AMA1.100.tif"), width = 8, height = 3, units = "in", res = 1200)
  par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))
  
  ggplot(PRISM.AMA1.100[1:48,], aes(x = rowid, y = X13_1.PfAMA1.100ug.ml_1)) + geom_point(col = "blue") + theme_bw() + 
    labs(x = "Row ID", y = "Normalized Log2(Positive/Negative)", title = "PRISM AMA1 100 µg/mL") + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 9)) +
    ylim(0,9)
  
  graphics.off()
  
#2. CP3 and PRISM without taking a ratio to the negative control pool 
  
  #prepare the final data frame without applying the subtractNeg function, yes Neg, no Swazi
  Finaldata <- rbind(CP3,PRISM, Neg)
  
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
  
  RawNorm.df <- rbind(PB1_meta, PB2_meta, PB3_meta)
  
  #need to change the sample_id to only CP3, PRISM, and Swazi 
  #so that the model knows they are the same sample
  RawNorm.df$sample <- "CP3"
  RawNorm.df$sample[c(grep("PRISM", RawNorm.df$sample_id))] <- "PRISM"
  RawNorm.df$sample[c(grep("Neg", RawNorm.df$sample_id))] <- "Neg"
  
  #the names of the columns have spaces. So they cannot be used in lmer and other functions
  newnames <- make.names(colnames(RawNorm.df))
  colnames(RawNorm.df) <- newnames
  
  #isolate data for CP3, PRISM, and Neg again
  N.CP3all <- filter(RawNorm.df, sample == "CP3")
  N.PRISMall <- filter(RawNorm.df, sample == "PRISM")
  N.Negall <- filter(RawNorm.df, sample == "Neg")

  N.PRISMall <- tibble::rowid_to_column(N.PRISMall)
  N.CP3all <- tibble::rowid_to_column(N.CP3all)
  N.Negall <- tibble::rowid_to_column(N.Negall)
  
  #sort/order by desired column, highest to lowest values
  N.CP3.AMA1.100 <- N.CP3all[order(N.CP3all$X13_1.PfAMA1.100ug.ml_1, decreasing = TRUE),]
  rowidorder <- c(as.character(N.CP3.AMA1.100$rowid))
  N.CP3.AMA1.100$rowid <- factor(N.CP3.AMA1.100$rowid, levels = unique(rowidorder))
  
  Neg.CP3.AMA1 <- N.Negall[order(N.CP3all$X13_1.PfAMA1.100ug.ml_1, decreasing = TRUE),]
  Neg.CP3.AMA1$rowid <- factor(Neg.CP3.AMA1$rowid, levels = unique(rowidorder))
  #bind Neg and CP3 to use as one data frame in plot
  Neg.CP3.AMA1 <- rbind(N.CP3.AMA1.100[1:48,], Neg.CP3.AMA1[1:48,])
  
  N.PRISM.AMA1.100 <- N.PRISMall[order(N.PRISMall$X13_1.PfAMA1.100ug.ml_1, decreasing = TRUE),]
  rowidorder <- c(as.character(N.PRISM.AMA1.100$rowid))
  N.PRISM.AMA1.100$rowid <- factor(N.PRISM.AMA1.100$rowid, levels = unique(rowidorder))
  
  Neg.PRISM.AMA1 <- N.Negall[order(N.PRISMall$X13_1.PfAMA1.100ug.ml_1, decreasing = TRUE),]
  Neg.PRISM.AMA1$rowid <- factor(Neg.PRISM.AMA1$rowid, levels = unique(rowidorder))
  #bind Neg and PRISM to use as one data frame in plot
  Neg.PRISM.AMA1 <- rbind(N.PRISM.AMA1.100[1:48,], Neg.PRISM.AMA1[1:48,])
  
  #find top 48 (10%) of conditions and export table of conditions matching top 10% rowid. 
  conditionsAMA1.100 <- merge(N.CP3.AMA1.100, conditions, sort = FALSE)
  conditionAMA1sub <- conditionsAMA1.100[1:48, 1:5]
  write.csv(conditionAMA1sub, file = "N.CP3top10AMA1.100.csv")
  
  N.PRISM.cond.AMA1.100 <- merge(N.PRISM.AMA1.100, conditions, sort = FALSE)
  N.PRISM.cond.AMA1sub <- N.PRISM.cond.AMA1.100[1:48, 1:5]
  write.csv(N.PRISM.cond.AMA1sub, file = "N.PRISMtop10AMA1.100.csv")
  
  #plot the top 10% values and rowids for CP3 alone and Neg
  png(filename = paste0("N.CP3.Top.AMA1.100.tif"), width = 8, height = 3, units = "in", res = 1200)
  par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))
  
  ggplot(Neg.CP3.AMA1, aes(x = rowid, y = X13_1.PfAMA1.100ug.ml_1, color = sample)) + geom_point() + theme_bw() + 
    labs(x = "Row ID", y = "Normalized Log2(MFI)", title = "CP3 AMA1 100 µg/mL") + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 9)) +
    ylim(0,10)
  
  graphics.off()
  
  #plot the top 10% values and rowids for PRISM alone and Negs
  #change factor levels of sample so that Neg is 2nd and PRISM is 1st
  Neg.PRISM.AMA1$sample <- factor(Neg.PRISM.AMA1$sample, levels = c("PRISM", "Neg"))
  
  png(filename = paste0("N.PRISM.Top.AMA1.100.tif"), width = 8, height = 3, units = "in", res = 1200)
  par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))
  
  ggplot(Neg.PRISM.AMA1, aes(x = rowid, y = X13_1.PfAMA1.100ug.ml_1, color = sample)) + geom_point() + theme_bw() + 
    labs(x = "Row ID", y = "Normalized Log2(MFI)", title = "PRISM AMA1 100 µg/mL") + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 9)) +
    ylim(0,10)
  
  graphics.off()
  
# MSP1.19 100 µg/mL  
  
  
  
  
  
  
 #old code for plotting positives and negatives on the same plots 
  png(filename = paste0("N.CP3.Top.AMA1.100.tif"), width = 8, height = 3, units = "in", res = 1200)
  par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))
  
  base <- ggplot(N.CP3.AMA1.100[1:48,], aes(x = rowid, y = X13_1.PfAMA1.100ug.ml_1, color = "CP3")) + geom_point(col = "red") + theme_bw() + 
    labs(x = "Row ID", y = "Normalized Log2(MFI)", title = "CP3 AMA1 100 µg/mL") + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 9)) +
    ylim(0,10) + scale_color_manual(name="Sample", values=c(CP3 ="red", Neg ="blue"))
  
  base + geom_point(data = Neg.CP3.AMA1[1:48,], aes(x=rowid, y=X13_1.PfAMA1.100ug.ml_1, color = "Neg"))
  
  graphics.off()
  
#old code for vertical plots
geomeanAMA1.100$rowid <- factor(geomeanAMA1.100$rowid, levels = rev(unique(rowidorder)))

png(filename = paste0("V.AMA1.100.GM.tif"), width = 3.5, height = 4.5, units = "in", res = 1200)
par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))

ggplot(geomeanAMA1.100[1:48,], aes(x = rowid, y = X13_1.PfAMA1.100ug.ml_1)) + theme_bw() + 
  geom_point() + theme(axis.text.y = element_text(size = 7)) +
  coord_flip() + labs(x = "Row ID", y = "Normalized Log2(Positive/Negative)", title = "AMA1 100 ?g/mL") +
  ylim(0,8)

graphics.off()

