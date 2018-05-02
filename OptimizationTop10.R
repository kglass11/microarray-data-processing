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
  
  #1. CP3 and PRISM with ratio from positive to negative (CP3all has everything for CP3 samples only)
  
  #sort/order by desired column, highest to lowest values
  colnames(CP3all[18])
  CP3.MSP1.19.100 <- CP3all[order(CP3all$X19_1.PfMSP1.19.100ug.ml_1, decreasing = TRUE),]
  rowidorder <- c(as.character(CP3.MSP1.19.100$rowid))
  CP3.MSP1.19.100$rowid <- factor(CP3.MSP1.19.100$rowid, levels = unique(rowidorder))
  
  PRISM.MSP1.19.100 <- PRISMall[order(PRISMall$X19_1.PfMSP1.19.100ug.ml_1, decreasing = TRUE),]
  rowidorder <- c(as.character(PRISM.MSP1.19.100$rowid))
  PRISM.MSP1.19.100$rowid <- factor(PRISM.MSP1.19.100$rowid, levels = unique(rowidorder))
  
  #find top 48 (10%) of conditions and export table of conditions matching top 10% rowid. 
  conditionsMSP1.19.100 <- merge(CP3.MSP1.19.100, conditions, sort = FALSE)
  conditionMSP1.19sub <- conditionsMSP1.19.100[1:48, 1:5]
  write.csv(conditionMSP1.19sub, file = "CP3top10MSP1.19.100.csv")
  
  PRISM.cond.MSP1.19.100 <- merge(PRISM.MSP1.19.100, conditions, sort = FALSE)
  PRISM.cond.MSP1.19sub <- PRISM.cond.MSP1.19.100[1:48, 1:5]
  write.csv(PRISM.cond.MSP1.19sub, file = "PRISMtop10MSP1.19.100.csv")
  
  #make one data frame with PRISM and CP3 best conditions
  BothMSP1.19.100 <- rbind(conditionsMSP1.19.100[1:48,], PRISM.cond.MSP1.19.100[1:48,])
  
  #plot the top 10% values and rowids for both CP3 and PRISM on one plot - this isn't working yet
  png(filename = paste0("Both.Top.MSP1.19.100.tif"), width = 8, height = 3, units = "in", res = 1200)
  par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))
  
  base <- ggplot(CP3.MSP1.19.100[1:48,]) + geom_point(aes(x = rowid, y = X19_1.PfMSP1.19.100ug.ml_1), col = "red") + theme_bw() + 
    labs(x = "Row ID", y = "Normalized Log2(Positive/Negative)", title = "Best Conditions MSP1.19 100 µg/mL") + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 9)) +
    ylim(0,9)
  
  base + geom_point(data = PRISM.MSP1.19.100[1:48,], aes(x = rowid, y = X13_1.PfAMA1.100ug.ml_1), col= "green")
  
  graphics.off()
  
  #plot the top 10% values and rowids for CP3 alone
  png(filename = paste0("CP3.Top.MSP1.19.100.tif"), width = 8, height = 3, units = "in", res = 1200)
  par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))
  
  ggplot(CP3.MSP1.19.100[1:48,], aes(x = rowid, y = X19_1.PfMSP1.19.100ug.ml_1)) + geom_point(col = "red") + theme_bw() + 
    labs(x = "Row ID", y = "Normalized Log2(Positive/Negative)", title = "CP3 MSP1.19 100 µg/mL") + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 9)) +
    ylim(0,9)
  
  graphics.off()
  
  #plot the top 10% values and rowids for PRISM alone
  png(filename = paste0("PRISM.Top.MSP1.19.100.tif"), width = 8, height = 3, units = "in", res = 1200)
  par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))
  
  ggplot(PRISM.MSP1.19.100[1:48,], aes(x = rowid, y = X19_1.PfMSP1.19.100ug.ml_1)) + geom_point(col = "blue") + theme_bw() + 
    labs(x = "Row ID", y = "Normalized Log2(Positive/Negative)", title = "PRISM MSP1.19 100 µg/mL") + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 9)) +
    ylim(0,9)
  
  graphics.off() 
  
 #2. without negative control ratio 
  
  #sort/order by desired column, highest to lowest values
  N.CP3.MSP1.19.100 <- N.CP3all[order(N.CP3all$X19_1.PfMSP1.19.100ug.ml_1, decreasing = TRUE),]
  rowidorder <- c(as.character(N.CP3.MSP1.19.100$rowid))
  N.CP3.MSP1.19.100$rowid <- factor(N.CP3.MSP1.19.100$rowid, levels = unique(rowidorder))
  
  Neg.CP3.MSP1.19 <- N.Negall[order(N.CP3all$X19_1.PfMSP1.19.100ug.ml_1, decreasing = TRUE),]
  Neg.CP3.MSP1.19$rowid <- factor(Neg.CP3.MSP1.19$rowid, levels = unique(rowidorder))
  #bind Neg and CP3 to use as one data frame in plot
  Neg.CP3.MSP1.19 <- rbind(N.CP3.MSP1.19.100[1:48,], Neg.CP3.MSP1.19[1:48,])
  
  N.PRISM.MSP1.19.100 <- N.PRISMall[order(N.PRISMall$X19_1.PfMSP1.19.100ug.ml_1, decreasing = TRUE),]
  rowidorder <- c(as.character(N.PRISM.MSP1.19.100$rowid))
  N.PRISM.MSP1.19.100$rowid <- factor(N.PRISM.MSP1.19.100$rowid, levels = unique(rowidorder))
  
  Neg.PRISM.MSP1.19 <- N.Negall[order(N.PRISMall$X19_1.PfMSP1.19.100ug.ml_1, decreasing = TRUE),]
  Neg.PRISM.MSP1.19$rowid <- factor(Neg.PRISM.MSP1.19$rowid, levels = unique(rowidorder))
  #bind Neg and PRISM to use as one data frame in plot
  Neg.PRISM.MSP1.19 <- rbind(N.PRISM.MSP1.19.100[1:48,], Neg.PRISM.MSP1.19[1:48,])
  
  #find top 48 (10%) of conditions and export table of conditions matching top 10% rowid. 
  conditionsMSP1.19.100 <- merge(N.CP3.MSP1.19.100, conditions, sort = FALSE)
  conditionMSP1.19sub <- conditionsMSP1.19.100[1:48, 1:5]
  write.csv(conditionMSP1.19sub, file = "N.CP3top10MSP1.19.100.csv")
  
  N.PRISM.cond.MSP1.19.100 <- merge(N.PRISM.MSP1.19.100, conditions, sort = FALSE)
  N.PRISM.cond.MSP1.19sub <- N.PRISM.cond.MSP1.19.100[1:48, 1:5]
  write.csv(N.PRISM.cond.MSP1.19sub, file = "N.PRISMtop10MSP1.19.100.csv")
  
  #plot the top 10% values and rowids for CP3 alone and Neg
  png(filename = paste0("N.CP3.Top.MSP1.19.100.tif"), width = 8, height = 3, units = "in", res = 1200)
  par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))
  
  ggplot(Neg.CP3.MSP1.19, aes(x = rowid, y = X19_1.PfMSP1.19.100ug.ml_1, color = sample)) + geom_point() + theme_bw() + 
    labs(x = "Row ID", y = "Normalized Log2(MFI)", title = "CP3 MSP1.19 100 µg/mL") + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 9)) +
    ylim(0,10)
  
  graphics.off()
  
  #plot the top 10% values and rowids for PRISM alone and Negs
  #change factor levels of sample so that Neg is 2nd and PRISM is 1st
  Neg.PRISM.MSP1.19$sample <- factor(Neg.PRISM.MSP1.19$sample, levels = c("PRISM", "Neg"))
  
  png(filename = paste0("N.PRISM.Top.MSP1.19.100.tif"), width = 8, height = 3, units = "in", res = 1200)
  par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))
  
  ggplot(Neg.PRISM.MSP1.19, aes(x = rowid, y = X19_1.PfMSP1.19.100ug.ml_1, color = sample)) + geom_point() + theme_bw() + 
    labs(x = "Row ID", y = "Normalized Log2(MFI)", title = "PRISM MSP1.19 100 µg/mL") + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 9)) +
    ylim(0,10)
  
  graphics.off()
  
  
# Hyp2 100 µg/mL  
  
  #1. CP3 and PRISM with ratio from positive to negative (CP3all has everything for CP3 samples only)
  
  #sort/order by desired column, highest to lowest values
  colnames(CP3all[24])
  CP3.Hyp2.100 <- CP3all[order(CP3all$X25_1.Hyp2.100ug.ml_1, decreasing = TRUE),]
  rowidorder <- c(as.character(CP3.Hyp2.100$rowid))
  CP3.Hyp2.100$rowid <- factor(CP3.Hyp2.100$rowid, levels = unique(rowidorder))
  
  PRISM.Hyp2.100 <- PRISMall[order(PRISMall$X25_1.Hyp2.100ug.ml_1, decreasing = TRUE),]
  rowidorder <- c(as.character(PRISM.Hyp2.100$rowid))
  PRISM.Hyp2.100$rowid <- factor(PRISM.Hyp2.100$rowid, levels = unique(rowidorder))
  
  #find top 48 (10%) of conditions and export table of conditions matching top 10% rowid. 
  conditionsHyp2.100 <- merge(CP3.Hyp2.100, conditions, sort = FALSE)
  conditionHyp2sub <- conditionsHyp2.100[1:48, 1:5]
  write.csv(conditionHyp2sub, file = "CP3top10Hyp2.100.csv")
  
  PRISM.cond.Hyp2.100 <- merge(PRISM.Hyp2.100, conditions, sort = FALSE)
  PRISM.cond.Hyp2sub <- PRISM.cond.Hyp2.100[1:48, 1:5]
  write.csv(PRISM.cond.Hyp2sub, file = "PRISMtop10Hyp2.100.csv")
  
  #make one data frame with PRISM and CP3 best conditions
  BothHyp2.100 <- rbind(conditionsHyp2.100[1:48,], PRISM.cond.Hyp2.100[1:48,])
  
  #plot the top 10% values and rowids for both CP3 and PRISM on one plot - this isn't working yet
  png(filename = paste0("Both.Top.Hyp2.100.tif"), width = 8, height = 3, units = "in", res = 1200)
  par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))
  
  base <- ggplot(CP3.Hyp2.100[1:48,]) + geom_point(aes(x = rowid, y = X25_1.Hyp2.100ug.ml_1), col = "red") + theme_bw() + 
    labs(x = "Row ID", y = "Normalized Log2(Positive/Negative)", title = "Best Conditions Hyp2 100 µg/mL") + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 9)) +
    ylim(0,9)
  
  base + geom_point(data = PRISM.Hyp2.100[1:48,], aes(x = rowid, y = X13_1.PfAMA1.100ug.ml_1), col= "green")
  
  graphics.off()
  
  #plot the top 10% values and rowids for CP3 alone
  png(filename = paste0("CP3.Top.Hyp2.100.tif"), width = 8, height = 3, units = "in", res = 1200)
  par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))
  
  ggplot(CP3.Hyp2.100[1:48,], aes(x = rowid, y = X25_1.Hyp2.100ug.ml_1)) + geom_point(col = "red") + theme_bw() + 
    labs(x = "Row ID", y = "Normalized Log2(Positive/Negative)", title = "CP3 Hyp2 100 µg/mL") + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 9)) +
    ylim(0,9)
  
  graphics.off()
  
  #plot the top 10% values and rowids for PRISM alone
  png(filename = paste0("PRISM.Top.Hyp2.100.tif"), width = 8, height = 3, units = "in", res = 1200)
  par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))
  
  ggplot(PRISM.Hyp2.100[1:48,], aes(x = rowid, y = X25_1.Hyp2.100ug.ml_1)) + geom_point(col = "blue") + theme_bw() + 
    labs(x = "Row ID", y = "Normalized Log2(Positive/Negative)", title = "PRISM Hyp2 100 µg/mL") + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 9)) +
    ylim(0,9)
  
  graphics.off() 
  
  #2. without negative control ratio 
  
  #sort/order by desired column, highest to lowest values
  N.CP3.Hyp2.100 <- N.CP3all[order(N.CP3all$X25_1.Hyp2.100ug.ml_1, decreasing = TRUE),]
  rowidorder <- c(as.character(N.CP3.Hyp2.100$rowid))
  N.CP3.Hyp2.100$rowid <- factor(N.CP3.Hyp2.100$rowid, levels = unique(rowidorder))
  
  Neg.CP3.Hyp2 <- N.Negall[order(N.CP3all$X25_1.Hyp2.100ug.ml_1, decreasing = TRUE),]
  Neg.CP3.Hyp2$rowid <- factor(Neg.CP3.Hyp2$rowid, levels = unique(rowidorder))
  #bind Neg and CP3 to use as one data frame in plot
  Neg.CP3.Hyp2 <- rbind(N.CP3.Hyp2.100[1:48,], Neg.CP3.Hyp2[1:48,])
  
  N.PRISM.Hyp2.100 <- N.PRISMall[order(N.PRISMall$X25_1.Hyp2.100ug.ml_1, decreasing = TRUE),]
  rowidorder <- c(as.character(N.PRISM.Hyp2.100$rowid))
  N.PRISM.Hyp2.100$rowid <- factor(N.PRISM.Hyp2.100$rowid, levels = unique(rowidorder))
  
  Neg.PRISM.Hyp2 <- N.Negall[order(N.PRISMall$X25_1.Hyp2.100ug.ml_1, decreasing = TRUE),]
  Neg.PRISM.Hyp2$rowid <- factor(Neg.PRISM.Hyp2$rowid, levels = unique(rowidorder))
  #bind Neg and PRISM to use as one data frame in plot
  Neg.PRISM.Hyp2 <- rbind(N.PRISM.Hyp2.100[1:48,], Neg.PRISM.Hyp2[1:48,])
  
  #find top 48 (10%) of conditions and export table of conditions matching top 10% rowid. 
  conditionsHyp2.100 <- merge(N.CP3.Hyp2.100, conditions, sort = FALSE)
  conditionHyp2sub <- conditionsHyp2.100[1:48, 1:5]
  write.csv(conditionHyp2sub, file = "N.CP3top10Hyp2.100.csv")
  
  N.PRISM.cond.Hyp2.100 <- merge(N.PRISM.Hyp2.100, conditions, sort = FALSE)
  N.PRISM.cond.Hyp2sub <- N.PRISM.cond.Hyp2.100[1:48, 1:5]
  write.csv(N.PRISM.cond.Hyp2sub, file = "N.PRISMtop10Hyp2.100.csv")
  
  #plot the top 10% values and rowids for CP3 alone and Neg
  png(filename = paste0("N.CP3.Top.Hyp2.100.tif"), width = 8, height = 3, units = "in", res = 1200)
  par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))
  
  ggplot(Neg.CP3.Hyp2, aes(x = rowid, y = X25_1.Hyp2.100ug.ml_1, color = sample)) + geom_point() + theme_bw() + 
    labs(x = "Row ID", y = "Normalized Log2(MFI)", title = "CP3 Hyp2 100 µg/mL") + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 9)) +
    ylim(0,10)
  
  graphics.off()
  
  #plot the top 10% values and rowids for PRISM alone and Negs
  #change factor levels of sample so that Neg is 2nd and PRISM is 1st
  Neg.PRISM.Hyp2$sample <- factor(Neg.PRISM.Hyp2$sample, levels = c("PRISM", "Neg"))
  
  png(filename = paste0("N.PRISM.Top.Hyp2.100.tif"), width = 8, height = 3, units = "in", res = 1200)
  par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))
  
  ggplot(Neg.PRISM.Hyp2, aes(x = rowid, y = X25_1.Hyp2.100ug.ml_1, color = sample)) + geom_point() + theme_bw() + 
    labs(x = "Row ID", y = "Normalized Log2(MFI)", title = "PRISM Hyp2 100 µg/mL") + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 9)) +
    ylim(0,10)
  
  graphics.off()
  
  
# EPF1v2 100 µg/mL  
  
  #1. CP3 and PRISM with ratio from positive to negative (CP3all has everything for CP3 samples only)
  
  #sort/order by desired column, highest to lowest values
  colnames(CP3all[36])
  CP3.EPF1v2.100 <- CP3all[order(CP3all$X37_1.EPF1v2.100ug.ml_1, decreasing = TRUE),]
  rowidorder <- c(as.character(CP3.EPF1v2.100$rowid))
  CP3.EPF1v2.100$rowid <- factor(CP3.EPF1v2.100$rowid, levels = unique(rowidorder))
  
  PRISM.EPF1v2.100 <- PRISMall[order(PRISMall$X37_1.EPF1v2.100ug.ml_1, decreasing = TRUE),]
  rowidorder <- c(as.character(PRISM.EPF1v2.100$rowid))
  PRISM.EPF1v2.100$rowid <- factor(PRISM.EPF1v2.100$rowid, levels = unique(rowidorder))
  
  #find top 48 (10%) of conditions and export table of conditions matching top 10% rowid. 
  conditionsEPF1v2.100 <- merge(CP3.EPF1v2.100, conditions, sort = FALSE)
  conditionEPF1v2sub <- conditionsEPF1v2.100[1:48, 1:5]
  write.csv(conditionEPF1v2sub, file = "CP3top10EPF1v2.100.csv")
  
  PRISM.cond.EPF1v2.100 <- merge(PRISM.EPF1v2.100, conditions, sort = FALSE)
  PRISM.cond.EPF1v2sub <- PRISM.cond.EPF1v2.100[1:48, 1:5]
  write.csv(PRISM.cond.EPF1v2sub, file = "PRISMtop10EPF1v2.100.csv")
  
  #make one data frame with PRISM and CP3 best conditions
  BothEPF1v2.100 <- rbind(conditionsEPF1v2.100[1:48,], PRISM.cond.EPF1v2.100[1:48,])
  
  #plot the top 10% values and rowids for both CP3 and PRISM on one plot - this isn't working yet
  png(filename = paste0("Both.Top.EPF1v2.100.tif"), width = 8, height = 3, units = "in", res = 1200)
  par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))
  
  base <- ggplot(CP3.EPF1v2.100[1:48,]) + geom_point(aes(x = rowid, y = X37_1.EPF1v2.100ug.ml_1), col = "red") + theme_bw() + 
    labs(x = "Row ID", y = "Normalized Log2(Positive/Negative)", title = "Best Conditions EPF1v2 100 µg/mL") + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 9)) +
    ylim(0,9)
  
  base + geom_point(data = PRISM.EPF1v2.100[1:48,], aes(x = rowid, y = X13_1.PfAMA1.100ug.ml_1), col= "green")
  
  graphics.off()
  
  #plot the top 10% values and rowids for CP3 alone
  png(filename = paste0("CP3.Top.EPF1v2.100.tif"), width = 8, height = 3, units = "in", res = 1200)
  par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))
  
  ggplot(CP3.EPF1v2.100[1:48,], aes(x = rowid, y = X37_1.EPF1v2.100ug.ml_1)) + geom_point(col = "red") + theme_bw() + 
    labs(x = "Row ID", y = "Normalized Log2(Positive/Negative)", title = "CP3 EPF1v2 100 µg/mL") + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 9)) +
    ylim(0,9)
  
  graphics.off()
  
  #plot the top 10% values and rowids for PRISM alone
  png(filename = paste0("PRISM.Top.EPF1v2.100.tif"), width = 8, height = 3, units = "in", res = 1200)
  par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))
  
  ggplot(PRISM.EPF1v2.100[1:48,], aes(x = rowid, y = X37_1.EPF1v2.100ug.ml_1)) + geom_point(col = "blue") + theme_bw() + 
    labs(x = "Row ID", y = "Normalized Log2(Positive/Negative)", title = "PRISM EPF1v2 100 µg/mL") + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 9)) +
    ylim(0,9)
  
  graphics.off() 
  
  #2. without negative control ratio 
  
  #sort/order by desired column, highest to lowest values
  N.CP3.EPF1v2.100 <- N.CP3all[order(N.CP3all$X37_1.EPF1v2.100ug.ml_1, decreasing = TRUE),]
  rowidorder <- c(as.character(N.CP3.EPF1v2.100$rowid))
  N.CP3.EPF1v2.100$rowid <- factor(N.CP3.EPF1v2.100$rowid, levels = unique(rowidorder))
  
  Neg.CP3.EPF1v2 <- N.Negall[order(N.CP3all$X37_1.EPF1v2.100ug.ml_1, decreasing = TRUE),]
  Neg.CP3.EPF1v2$rowid <- factor(Neg.CP3.EPF1v2$rowid, levels = unique(rowidorder))
  #bind Neg and CP3 to use as one data frame in plot
  Neg.CP3.EPF1v2 <- rbind(N.CP3.EPF1v2.100[1:48,], Neg.CP3.EPF1v2[1:48,])
  
  N.PRISM.EPF1v2.100 <- N.PRISMall[order(N.PRISMall$X37_1.EPF1v2.100ug.ml_1, decreasing = TRUE),]
  rowidorder <- c(as.character(N.PRISM.EPF1v2.100$rowid))
  N.PRISM.EPF1v2.100$rowid <- factor(N.PRISM.EPF1v2.100$rowid, levels = unique(rowidorder))
  
  Neg.PRISM.EPF1v2 <- N.Negall[order(N.PRISMall$X37_1.EPF1v2.100ug.ml_1, decreasing = TRUE),]
  Neg.PRISM.EPF1v2$rowid <- factor(Neg.PRISM.EPF1v2$rowid, levels = unique(rowidorder))
  #bind Neg and PRISM to use as one data frame in plot
  Neg.PRISM.EPF1v2 <- rbind(N.PRISM.EPF1v2.100[1:48,], Neg.PRISM.EPF1v2[1:48,])
  
  #find top 48 (10%) of conditions and export table of conditions matching top 10% rowid. 
  conditionsEPF1v2.100 <- merge(N.CP3.EPF1v2.100, conditions, sort = FALSE)
  conditionEPF1v2sub <- conditionsEPF1v2.100[1:48, 1:5]
  write.csv(conditionEPF1v2sub, file = "N.CP3top10EPF1v2.100.csv")
  
  N.PRISM.cond.EPF1v2.100 <- merge(N.PRISM.EPF1v2.100, conditions, sort = FALSE)
  N.PRISM.cond.EPF1v2sub <- N.PRISM.cond.EPF1v2.100[1:48, 1:5]
  write.csv(N.PRISM.cond.EPF1v2sub, file = "N.PRISMtop10EPF1v2.100.csv")
  
  #plot the top 10% values and rowids for CP3 alone and Neg
  png(filename = paste0("N.CP3.Top.EPF1v2.100.tif"), width = 8, height = 3, units = "in", res = 1200)
  par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))
  
  ggplot(Neg.CP3.EPF1v2, aes(x = rowid, y = X37_1.EPF1v2.100ug.ml_1, color = sample)) + geom_point() + theme_bw() + 
    labs(x = "Row ID", y = "Normalized Log2(MFI)", title = "CP3 EPF1v2 100 µg/mL") + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 9)) +
    ylim(0,10)
  
  graphics.off()
  
  #plot the top 10% values and rowids for PRISM alone and Negs
  #change factor levels of sample so that Neg is 2nd and PRISM is 1st
  Neg.PRISM.EPF1v2$sample <- factor(Neg.PRISM.EPF1v2$sample, levels = c("PRISM", "Neg"))
  
  png(filename = paste0("N.PRISM.Top.EPF1v2.100.tif"), width = 8, height = 3, units = "in", res = 1200)
  par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))
  
  ggplot(Neg.PRISM.EPF1v2, aes(x = rowid, y = X37_1.EPF1v2.100ug.ml_1, color = sample)) + geom_point() + theme_bw() + 
    labs(x = "Row ID", y = "Normalized Log2(MFI)", title = "PRISM EPF1v2 100 µg/mL") + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 9)) +
    ylim(0,10)
  
  graphics.off()
  
# GexP18 100 �g/mL  
  
  #1. CP3 and PRISM with ratio from positive to negative (CP3all has everything for CP3 samples only)
  
  #sort/order by desired column, highest to lowest values
  colnames(CP3all[30])
  CP3.GEXP18.100 <- CP3all[order(CP3all$X31_1.GEXP18.100ug.ml_1, decreasing = TRUE),]
  rowidorder <- c(as.character(CP3.GEXP18.100$rowid))
  CP3.GEXP18.100$rowid <- factor(CP3.GEXP18.100$rowid, levels = unique(rowidorder))
  
  PRISM.GEXP18.100 <- PRISMall[order(PRISMall$X31_1.GEXP18.100ug.ml_1, decreasing = TRUE),]
  rowidorder <- c(as.character(PRISM.GEXP18.100$rowid))
  PRISM.GEXP18.100$rowid <- factor(PRISM.GEXP18.100$rowid, levels = unique(rowidorder))
  
  #find top 48 (10%) of conditions and export table of conditions matching top 10% rowid. 
  conditionsGEXP18.100 <- merge(CP3.GEXP18.100, conditions, sort = FALSE)
  conditionGEXP18sub <- conditionsGEXP18.100[1:48, 1:5]
  write.csv(conditionGEXP18sub, file = "CP3top10GEXP18.100.csv")
  
  PRISM.cond.GEXP18.100 <- merge(PRISM.GEXP18.100, conditions, sort = FALSE)
  PRISM.cond.GEXP18sub <- PRISM.cond.GEXP18.100[1:48, 1:5]
  write.csv(PRISM.cond.GEXP18sub, file = "PRISMtop10GEXP18.100.csv")
  
  #make one data frame with PRISM and CP3 best conditions
  BothGEXP18.100 <- rbind(conditionsGEXP18.100[1:48,], PRISM.cond.GEXP18.100[1:48,])
  
  #plot the top 10% values and rowids for both CP3 and PRISM on one plot - this isn't working yet
  png(filename = paste0("Both.Top.GEXP18.100.tif"), width = 8, height = 3, units = "in", res = 1200)
  par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))
  
  base <- ggplot(CP3.GEXP18.100[1:48,]) + geom_point(aes(x = rowid, y = X31_1.GEXP18.100ug.ml_1), col = "red") + theme_bw() + 
    labs(x = "Row ID", y = "Normalized Log2(Positive/Negative)", title = "Best Conditions GEXP18 100 �g/mL") + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 9)) +
    ylim(0,9)
  
  base + geom_point(data = PRISM.GEXP18.100[1:48,], aes(x = rowid, y = X13_1.PfAMA1.100ug.ml_1), col= "green")
  
  graphics.off()
  
  #plot the top 10% values and rowids for CP3 alone
  png(filename = paste0("CP3.Top.GEXP18.100.tif"), width = 8, height = 3, units = "in", res = 1200)
  par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))
  
  ggplot(CP3.GEXP18.100[1:48,], aes(x = rowid, y = X31_1.GEXP18.100ug.ml_1)) + geom_point(col = "red") + theme_bw() + 
    labs(x = "Row ID", y = "Normalized Log2(Positive/Negative)", title = "CP3 GEXP18 100 �g/mL") + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 9)) +
    ylim(0,9)
  
  graphics.off()
  
  #plot the top 10% values and rowids for PRISM alone
  png(filename = paste0("PRISM.Top.GEXP18.100.tif"), width = 8, height = 3, units = "in", res = 1200)
  par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))
  
  ggplot(PRISM.GEXP18.100[1:48,], aes(x = rowid, y = X31_1.GEXP18.100ug.ml_1)) + geom_point(col = "blue") + theme_bw() + 
    labs(x = "Row ID", y = "Normalized Log2(Positive/Negative)", title = "PRISM GEXP18 100 �g/mL") + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 9)) +
    ylim(0,9)
  
  graphics.off() 
  
  
#GST 100 ug/mL 
  
  #1. CP3 and PRISM with ratio from positive to negative (CP3all has everything for CP3 samples only)
  
  #sort/order by desired column, highest to lowest values
  colnames(CP3all[42])
  CP3.GST.100 <- CP3all[order(CP3all$X43_1.GST.100ug.ml_1, decreasing = TRUE),]
  rowidorder <- c(as.character(CP3.GST.100$rowid))
  CP3.GST.100$rowid <- factor(CP3.GST.100$rowid, levels = unique(rowidorder))
  
  PRISM.GST.100 <- PRISMall[order(PRISMall$X43_1.GST.100ug.ml_1, decreasing = TRUE),]
  rowidorder <- c(as.character(PRISM.GST.100$rowid))
  PRISM.GST.100$rowid <- factor(PRISM.GST.100$rowid, levels = unique(rowidorder))
  
  #find top 48 (10%) of conditions and export table of conditions matching top 10% rowid. 
  conditionsGST.100 <- merge(CP3.GST.100, conditions, sort = FALSE)
  conditionGSTsub <- conditionsGST.100[1:48, 1:5]
  write.csv(conditionGSTsub, file = "CP3top10GST.100.csv")
  
  PRISM.cond.GST.100 <- merge(PRISM.GST.100, conditions, sort = FALSE)
  PRISM.cond.GSTsub <- PRISM.cond.GST.100[1:48, 1:5]
  write.csv(PRISM.cond.GSTsub, file = "PRISMtop10GST.100.csv")
  
  #make one data frame with PRISM and CP3 best conditions
  BothGST.100 <- rbind(conditionsGST.100[1:48,], PRISM.cond.GST.100[1:48,])
  
  #plot the top 10% values and rowids for both CP3 and PRISM on one plot - this isn't working yet
  png(filename = paste0("Both.Top.GST.100.tif"), width = 8, height = 3, units = "in", res = 1200)
  par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))
  
  base <- ggplot(CP3.GST.100[1:48,]) + geom_point(aes(x = rowid, y = X43_1.GST.100ug.ml_1), col = "red") + theme_bw() + 
    labs(x = "Row ID", y = "Normalized Log2(Positive/Negative)", title = "Best Conditions GST 100 �g/mL") + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 9)) +
    ylim(0,9)
  
  base + geom_point(data = PRISM.GST.100[1:48,], aes(x = rowid, y = X13_1.PfAMA1.100ug.ml_1), col= "green")
  
  graphics.off()
  
  #plot the top 10% values and rowids for CP3 alone
  png(filename = paste0("CP3.Top.GST.100.tif"), width = 8, height = 3, units = "in", res = 1200)
  par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))
  
  ggplot(CP3.GST.100[1:48,], aes(x = rowid, y = X43_1.GST.100ug.ml_1)) + geom_point(col = "red") + theme_bw() + 
    labs(x = "Row ID", y = "Normalized Log2(Positive/Negative)", title = "CP3 GST 100 �g/mL") + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 9)) +
    ylim(0,9)
  
  graphics.off()
  
  #plot the top 10% values and rowids for PRISM alone
  png(filename = paste0("PRISM.Top.GST.100.tif"), width = 8, height = 3, units = "in", res = 1200)
  par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))
  
  ggplot(PRISM.GST.100[1:48,], aes(x = rowid, y = X43_1.GST.100ug.ml_1)) + geom_point(col = "blue") + theme_bw() + 
    labs(x = "Row ID", y = "Normalized Log2(Positive/Negative)", title = "PRISM GST 100 �g/mL") + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 9)) +
    ylim(0,9)
  
  graphics.off() 
  
#Plot CP3 only, for all antigens on the same plot, all ranked highest to lowest pos/neg signal
  #x axis will be 1-48 for all, but this will no longer represent specific conditions
  
  #make a data frame with all antigens, 100 ug/mL, top 48 conditions, data for that antigen only
  #drop the rowid, because have that stored elsewhere
  top10.all.100 <- as.data.frame(cbind(CP3.AMA1.100[1:48,12],CP3.MSP1.19.100[1:48,18],CP3.Hyp2.100[1:48,24], CP3.GEXP18.100[1:48,30], CP3.EPF1v2.100[1:48,36],CP3.GST.100[1:48,42]))
  colnames(top10.all.100) <- as.character(c("AMA1", "MSP1-19", "Hyp2", "GEXP18", "EPF1v2", "GST"))
  
  #add another column with 1-48
  top10.all.100 <- tibble::rowid_to_column(top10.all.100)
  top10.all.100$rowid <- as.factor(top10.all.100$rowid)
  
  #melt the data
  top10.all.melt <- melt(top10.all.100)
  colnames(top10.all.melt) <- as.character(c("Number", "Antigen", "Value"))
  
  #set factor order to 1-48
  rowidorder <- c(as.character(top10.all.melt$rowid))
  top10.all.melt$rowid <- factor(top10.all.melt$rowid, levels = unique(rowidorder))
  
  #plot with color = antigen
  png(filename = paste0("CP3.Top10.ALL.100.tif"), width = 8, height = 3, units = "in", res = 1200)
  par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))
  
  ggplot(top10.all.melt, aes(x = Number, y = Value, color=Antigen)) + geom_point(shape=23) + theme_bw() + 
    labs(y = "Normalized Log2(Positive/Negative)", title = "CP3 values for top 10% of conditions for all antigens, 100 �g/mL") + 
    theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 8)) +
    ylim(0,9)
  
  graphics.off()
  
  #take out GST
  top10.noGST <- top10.all.melt[!top10.all.melt$Antigen == "GST",]
  antigenorder <- sort(c(as.character((top10.noGST$Antigen))))
  top10.noGST$Antigen <- factor(top10.noGST$Antigen, levels = unique(antigenorder))
  
  #plot without GST!!!
  png(filename = paste0("CP3.Top10.noGST.100.tif"), width = 8, height = 3, units = "in", res = 1200)
  par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))
  
  ggplot(top10.noGST, aes(x = Number, y = Value, color=Antigen)) + geom_point(shape=18) + theme_bw() + 
    labs(y = "Normalized Log2(Positive/Negative)", title = "CP3 values for top 10% of conditions, 100 �g/mL") + 
    theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 8)) +
    ylim(0,9)
  
  graphics.off()
  
  ###Find the conditions which are present for all antigens!! 
  #Ignoring the case without the negative ratio because the plots look terrible
  
  #First find the ones which are present in both CP3 and PRISM for each antigen
  bothAMA1 <- merge(conditionAMA1sub, PRISM.cond.AMA1sub, sort = FALSE) # one condition
  bothMSP1.19 <- merge(conditionMSP1.19sub, PRISM.cond.MSP1.19sub, sort = FALSE) #9 conditions
  bothHyp2 <- merge(conditionHyp2sub, PRISM.cond.Hyp2sub, sort = FALSE) #14 conditions
  bothEPF1v2 <- merge(conditionEPF1v2sub, PRISM.cond.EPF1v2sub, sort = FALSE) #12 conditions
  
  #check if bothAMA1 (only 1 condition :( ) is present in all antigens
  ALL <- merge(bothAMA1, bothMSP1.19) #NOT in MSP1.19
  ALL <- merge(bothAMA1, bothHyp2) #YES in Hyp2
  ALL <- merge(bothAMA1, bothEPF1v2) #NOT in EPF1v2
  
  #check if other 3 antigens have common top conditions
  ALL <- merge(bothMSP1.19, bothHyp2) #6 conditions in common
  ALL.3 <- merge(ALL, bothEPF1v2) #no conditions in common
  ALL <- merge(bothMSP1.19, bothEPF1v2) #1 condition in common
  ALL <- merge(bothHyp2, bothEPF1v2) #3 conditions in common
  
#check CP3 alone for common top conditions among all antigens 
  CP3alone2 <- merge(conditionAMA1sub, conditionMSP1.19sub) #28 conditions
  CP3alone3 <- merge(conditionHyp2sub, CP3alone2) #15 conditions
  CP3alone4 <- merge(conditionEPF1v2sub, CP3alone3) #7 conditions
  write.csv(CP3alone4, file = "CP3Top10ALLantigens.csv")

#check PRISM alone for common top conditions among all antigens
  PRISMalone2 <- merge(PRISM.cond.AMA1sub, PRISM.cond.MSP1.19sub) #15 conditions
  PRISMalone3 <- merge(PRISMalone2, PRISM.cond.Hyp2sub) #2 conditions
  PRISMalone4 <- merge(PRISMalone3, PRISM.cond.EPF1v2sub) #1 condition
  write.csv(PRISMalone4, file = "PRISMTop10ALLantigens.csv")
  
#check that PRISM and CP3 condition is not mutual
  All <- merge(PRISMalone4, CP3alone4) #0 conditions
  
#Check if conditions which came up for PRISM and CP3 alone are included in geomean selection
  All <- merge(CP3alone4, four) #0 conditions
  All <- merge(PRISMalone4, four) #0 conditions
  
   
#  #old code for plotting positives and negatives on the same plots 
#   png(filename = paste0("N.CP3.Top.AMA1.100.tif"), width = 8, height = 3, units = "in", res = 1200)
#   par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))
#   
#   base <- ggplot(N.CP3.AMA1.100[1:48,], aes(x = rowid, y = X13_1.PfAMA1.100ug.ml_1, color = "CP3")) + geom_point(col = "red") + theme_bw() + 
#     labs(x = "Row ID", y = "Normalized Log2(MFI)", title = "CP3 AMA1 100 µg/mL") + 
#     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 9)) +
#     ylim(0,10) + scale_color_manual(name="Sample", values=c(CP3 ="red", Neg ="blue"))
#   
#   base + geom_point(data = Neg.CP3.AMA1[1:48,], aes(x=rowid, y=X13_1.PfAMA1.100ug.ml_1, color = "Neg"))
#   
#   graphics.off()
#   
# #old code for vertical plots
# geomeanAMA1.100$rowid <- factor(geomeanAMA1.100$rowid, levels = rev(unique(rowidorder)))
# 
# png(filename = paste0("V.AMA1.100.GM.tif"), width = 3.5, height = 4.5, units = "in", res = 1200)
# par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))
# 
# ggplot(geomeanAMA1.100[1:48,], aes(x = rowid, y = X13_1.PfAMA1.100ug.ml_1)) + theme_bw() + 
#   geom_point() + theme(axis.text.y = element_text(size = 7)) +
#   coord_flip() + labs(x = "Row ID", y = "Normalized Log2(Positive/Negative)", title = "AMA1 100 ?g/mL") +
#   ylim(0,8)
# 
# graphics.off()
# 
