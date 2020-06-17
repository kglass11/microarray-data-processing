#continuing on from optimizationanalysis.R

load("OptimizationLMMready.RData")

#Updated 8/5/18
#This is an updated version which only has the analysis for the positive/negative ratio method

#Top 10% conditions based on signal from CP3 and PRISM separately,
  #all have GST subtraction the same way
  #all ratio of positive / negative

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
    labs(x = "Row ID", y = "Normalized Log2(Positive/Negative)", title = "Best Conditions AMA1 100 ????g/mL") + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 9)) +
    ylim(0,9)
  
  base + geom_point(data = PRISM.AMA1.100[1:48,], aes(x = rowid, y = X13_1.PfAMA1.100ug.ml_1), col= "green")
  
  graphics.off()
  
  #plot the top 10% values and rowids for CP3 alone
  png(filename = paste0("CP3.Top.AMA1.100.tif"), width = 8, height = 3, units = "in", res = 1200)
  par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))

  ggplot(CP3.AMA1.100[1:48,], aes(x = rowid, y = X13_1.PfAMA1.100ug.ml_1)) + geom_point(col = "red") + theme_bw() + 
  labs(x = "Row ID", y = "Normalized Log2(Positive/Negative)", title = "CP3 AMA1 100 ????g/mL") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 9)) +
  ylim(0,9)

  graphics.off()
  
  #plot the top 10% values and rowids for PRISM alone
  png(filename = paste0("PRISM.Top.AMA1.100.tif"), width = 8, height = 3, units = "in", res = 1200)
  par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))
  
  ggplot(PRISM.AMA1.100[1:48,], aes(x = rowid, y = X13_1.PfAMA1.100ug.ml_1)) + geom_point(col = "blue") + theme_bw() + 
    labs(x = "Row ID", y = "Normalized Log2(Positive/Negative)", title = "PRISM AMA1 100 ????g/mL") + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 9)) +
    ylim(0,9)
  
  graphics.off()
  
# MSP1.19 100 ????g/mL  
  
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
    labs(x = "Row ID", y = "Normalized Log2(Positive/Negative)", title = "Best Conditions MSP1.19 100 ????g/mL") + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 9)) +
    ylim(0,9)
  
  base + geom_point(data = PRISM.MSP1.19.100[1:48,], aes(x = rowid, y = X13_1.PfAMA1.100ug.ml_1), col= "green")
  
  graphics.off()
  
  #plot the top 10% values and rowids for CP3 alone
  png(filename = paste0("CP3.Top.MSP1.19.100.tif"), width = 8, height = 3, units = "in", res = 1200)
  par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))
  
  ggplot(CP3.MSP1.19.100[1:48,], aes(x = rowid, y = X19_1.PfMSP1.19.100ug.ml_1)) + geom_point(col = "red") + theme_bw() + 
    labs(x = "Row ID", y = "Normalized Log2(Positive/Negative)", title = "CP3 MSP1.19 100 ????g/mL") + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 9)) +
    ylim(0,9)
  
  graphics.off()
  
  #plot the top 10% values and rowids for PRISM alone
  png(filename = paste0("PRISM.Top.MSP1.19.100.tif"), width = 8, height = 3, units = "in", res = 1200)
  par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))
  
  ggplot(PRISM.MSP1.19.100[1:48,], aes(x = rowid, y = X19_1.PfMSP1.19.100ug.ml_1)) + geom_point(col = "blue") + theme_bw() + 
    labs(x = "Row ID", y = "Normalized Log2(Positive/Negative)", title = "PRISM MSP1.19 100 ????g/mL") + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 9)) +
    ylim(0,9)
  
  graphics.off() 
  
# Hyp2 100 ????g/mL  
  
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
    labs(x = "Row ID", y = "Normalized Log2(Positive/Negative)", title = "Best Conditions Hyp2 100 ????g/mL") + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 9)) +
    ylim(0,9)
  
  base + geom_point(data = PRISM.Hyp2.100[1:48,], aes(x = rowid, y = X13_1.PfAMA1.100ug.ml_1), col= "green")
  
  graphics.off()
  
  #plot the top 10% values and rowids for CP3 alone
  png(filename = paste0("CP3.Top.Hyp2.100.tif"), width = 8, height = 3, units = "in", res = 1200)
  par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))
  
  ggplot(CP3.Hyp2.100[1:48,], aes(x = rowid, y = X25_1.Hyp2.100ug.ml_1)) + geom_point(col = "red") + theme_bw() + 
    labs(x = "Row ID", y = "Normalized Log2(Positive/Negative)", title = "CP3 Hyp2 100 ????g/mL") + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 9)) +
    ylim(0,9)
  
  graphics.off()
  
  #plot the top 10% values and rowids for PRISM alone
  png(filename = paste0("PRISM.Top.Hyp2.100.tif"), width = 8, height = 3, units = "in", res = 1200)
  par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))
  
  ggplot(PRISM.Hyp2.100[1:48,], aes(x = rowid, y = X25_1.Hyp2.100ug.ml_1)) + geom_point(col = "blue") + theme_bw() + 
    labs(x = "Row ID", y = "Normalized Log2(Positive/Negative)", title = "PRISM Hyp2 100 ????g/mL") + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 9)) +
    ylim(0,9)
  
  graphics.off() 
  
# EPF1v2 100 ????g/mL  
  
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
    labs(x = "Row ID", y = "Normalized Log2(Positive/Negative)", title = "Best Conditions EPF1v2 100 ????g/mL") + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 9)) +
    ylim(0,9)
  
  base + geom_point(data = PRISM.EPF1v2.100[1:48,], aes(x = rowid, y = X13_1.PfAMA1.100ug.ml_1), col= "green")
  
  graphics.off()
  
  #plot the top 10% values and rowids for CP3 alone
  png(filename = paste0("CP3.Top.EPF1v2.100.tif"), width = 8, height = 3, units = "in", res = 1200)
  par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))
  
  ggplot(CP3.EPF1v2.100[1:48,], aes(x = rowid, y = X37_1.EPF1v2.100ug.ml_1)) + geom_point(col = "red") + theme_bw() + 
    labs(x = "Row ID", y = "Normalized Log2(Positive/Negative)", title = "CP3 EPF1v2 100 ????g/mL") + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 9)) +
    ylim(0,9)
  
  graphics.off()
  
  #plot the top 10% values and rowids for PRISM alone
  png(filename = paste0("PRISM.Top.EPF1v2.100.tif"), width = 8, height = 3, units = "in", res = 1200)
  par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))
  
  ggplot(PRISM.EPF1v2.100[1:48,], aes(x = rowid, y = X37_1.EPF1v2.100ug.ml_1)) + geom_point(col = "blue") + theme_bw() + 
    labs(x = "Row ID", y = "Normalized Log2(Positive/Negative)", title = "PRISM EPF1v2 100 ????g/mL") + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 9)) +
    ylim(0,9)
  
  graphics.off() 
  
  
# GexP18 100 ????g/mL  
  
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
    labs(x = "Row ID", y = "Normalized Log2(Positive/Negative)", title = "Best Conditions GEXP18 100 ?g/mL") + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 9)) +
    ylim(0,9)
  
  base + geom_point(data = PRISM.GEXP18.100[1:48,], aes(x = rowid, y = X13_1.PfAMA1.100ug.ml_1), col= "green")
  
  graphics.off()
  
  #plot the top 10% values and rowids for CP3 alone
  png(filename = paste0("CP3.Top.GEXP18.100.tif"), width = 8, height = 3, units = "in", res = 1200)
  par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))
  
  ggplot(CP3.GEXP18.100[1:48,], aes(x = rowid, y = X31_1.GEXP18.100ug.ml_1)) + geom_point(col = "red") + theme_bw() + 
    labs(x = "Row ID", y = "Normalized Log2(Positive/Negative)", title = "CP3 GEXP18 100 ?g/mL") + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 9)) +
    ylim(0,9)
  
  graphics.off()
  
  #plot the top 10% values and rowids for PRISM alone
  png(filename = paste0("PRISM.Top.GEXP18.100.tif"), width = 8, height = 3, units = "in", res = 1200)
  par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))
  
  ggplot(PRISM.GEXP18.100[1:48,], aes(x = rowid, y = X31_1.GEXP18.100ug.ml_1)) + geom_point(col = "blue") + theme_bw() + 
    labs(x = "Row ID", y = "Normalized Log2(Positive/Negative)", title = "PRISM GEXP18 100 ?g/mL") + 
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
    labs(x = "Row ID", y = "Normalized Log2(Positive/Negative)", title = "Best Conditions GST 100 ?g/mL") + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 9)) +
    ylim(0,9)
  
  base + geom_point(data = PRISM.GST.100[1:48,], aes(x = rowid, y = X13_1.PfAMA1.100ug.ml_1), col= "green")
  
  graphics.off()
  
  #plot the top 10% values and rowids for CP3 alone
  png(filename = paste0("CP3.Top.GST.100.tif"), width = 8, height = 3, units = "in", res = 1200)
  par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))
  
  ggplot(CP3.GST.100[1:48,], aes(x = rowid, y = X43_1.GST.100ug.ml_1)) + geom_point(col = "red") + theme_bw() + 
    labs(x = "Row ID", y = "Normalized Log2(Positive/Negative)", title = "CP3 GST 100 ?g/mL") + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 9)) +
    ylim(0,9)
  
  graphics.off()
  
  #plot the top 10% values and rowids for PRISM alone
  png(filename = paste0("PRISM.Top.GST.100.tif"), width = 8, height = 3, units = "in", res = 1200)
  par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))
  
  ggplot(PRISM.GST.100[1:48,], aes(x = rowid, y = X43_1.GST.100ug.ml_1)) + geom_point(col = "blue") + theme_bw() + 
    labs(x = "Row ID", y = "Normalized Log2(Positive/Negative)", title = "PRISM GST 100 ?g/mL") + 
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
  
  #plot with color = antigen
  png(filename = paste0("CP3.Top10.ALL.100.tif"), width = 8, height = 3, units = "in", res = 1200)
  par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))
  
  ggplot(top10.all.melt, aes(x = Number, y = Value, color=Antigen)) + geom_point(shape=23) + theme_bw() + 
    labs(y = "Normalized Log2(Positive/Negative)", title = "CP3 values for top 10% of conditions for all antigens, 100 ?g/mL") + 
    theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 8)) +
    ylim(0,9)
  
  graphics.off()
  
  #take out GST
  top10.noGST <- top10.all.melt[!top10.all.melt$Antigen == "GST",]
  antigenorder <- sort(c(as.character((top10.noGST$Antigen))))
  top10.noGST$Antigen <- factor(top10.noGST$Antigen, levels = unique(antigenorder))
  
  #plot without GST!!!
  png(filename = paste0("CP3.Top10.noGST.100.tif"), width = 8, height = 3.3, units = "in", res = 1200)
  
  ggplot(top10.noGST, aes(x = Number, y = Value, color=Antigen)) + geom_point(shape = 18, size = 2) + theme_bw() + 
    labs(y = "Normalized Log2(Positive/Negative)") + 
    theme(axis.text.x = element_text(angle = 90, size = 8.5, vjust = 0.5, color = "black")) +
    theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank()) +
    ylim(0,9)
  
  graphics.off()
  
  ###Find the conditions which are present for all antigens!! 
  #Ignoring the case without the negative ratio because the plots look terrible
  
  #First find the ones which are present in both CP3 and PRISM for each antigen
  bothAMA1 <- merge(conditionAMA1sub, PRISM.cond.AMA1sub, sort = FALSE) # 14 conditions
  bothMSP1.19 <- merge(conditionMSP1.19sub, PRISM.cond.MSP1.19sub, sort = FALSE) #14 conditions
  bothHyp2 <- merge(conditionHyp2sub, PRISM.cond.Hyp2sub, sort = FALSE) #16 conditions
  bothEPF1v2 <- merge(conditionEPF1v2sub, PRISM.cond.EPF1v2sub, sort = FALSE) #19 conditions
  
  #check if any conditions common to both PRISM and CP3 and 4 antigens
  ALL.2 <- merge(bothAMA1, bothMSP1.19) #7 in common
  ALL.3 <- merge(ALL.2, bothHyp2) #none in common
  ALL.3.E <- merge(ALL.2, bothEPF1v2) #none in common
  
  #check if other 3 antigens have common top conditions (without AMA1)
  ALL <- merge(bothMSP1.19, bothHyp2) #3 conditions in common
  ALL.3.3 <- merge(ALL, bothEPF1v2) #3 conditions in common
  
#check CP3 alone for common top conditions among all antigens 
  CP3alone2 <- merge(conditionAMA1sub, conditionMSP1.19sub) #28 conditions
  CP3alone3 <- merge(conditionHyp2sub, CP3alone2) #14 conditions
  CP3alone4 <- merge(conditionEPF1v2sub, CP3alone3) #8 conditions
  write.csv(CP3alone4, file = "CP3Top10ALLantigens.csv")

#check PRISM alone for common top conditions among all antigens
  PRISMalone2 <- merge(PRISM.cond.AMA1sub, PRISM.cond.MSP1.19sub) #15 conditions
  PRISMalone3 <- merge(PRISMalone2, PRISM.cond.Hyp2sub) #2 conditions
  PRISMalone4 <- merge(PRISMalone3, PRISM.cond.EPF1v2sub) #1 condition
  write.csv(PRISMalone4, file = "PRISMTop10ALLantigens.csv")
  
#check that PRISM and CP3 condition is not mutual
  All <- merge(PRISMalone4, CP3alone4) #0 conditions
  
#Check if conditions which came up for PRISM and CP3 alone are included in geomean selection
  All.CP3 <- merge(CP3alone4, four) #6 conditions
  All.PRISM <- merge(PRISMalone4, four) #0 conditions
  write.csv(All.CP3, file = "CP3andGM.ALLantigens.csv")
  
#Plot CP3 data for all antigens for the 8 common conditions ("best" conditions)
  
  #merge the conditions with the data to plot
  hello <- merge(CP3alone4, CP3all, sort = TRUE)
  
  #melt and subset data for the plot
  hellomelt <- melt(hello, variable.name = "Antigen")
  hellosub <- filter(hellomelt, Antigen == "X37_1.EPF1v2.100ug.ml_1" | Antigen == "X25_1.Hyp2.100ug.ml_1"
                     | Antigen == "X19_1.PfMSP1.19.100ug.ml_1" | Antigen == "X13_1.PfAMA1.100ug.ml_1")
  
  #vertical plot of CP3 for "best" 8 conditions - organized by highest to lowest minimum value
  png(filename = paste0("V.CP3.Common8.tif"), width = 4.2, height = 2.7, units = "in", res = 1200)
  
  ggplot(hellosub, aes(x = reorder(rowid, value, min), y = value, color = Antigen)) + theme_bw() + 
    geom_point(shape = 18, size = 2) + theme(axis.text.y = element_text(size = 10)) +
    scale_color_hue(labels = c("AMA1", "MSP1-19", "Hyp2", "EPF1v2")) +
    labs(x = "Row ID", y = "Normalized Log2(Positive/Negative)") +
    theme(axis.text.x = element_text(color = "black"), panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank()) +
    ylim(0,9)  + coord_flip()
  
  graphics.off()
  
  #Need to save the factor order for these row IDs so that can use the same factor order for all plots of top 8
  #not sure how to do the factor order outside the plot...reorder function was not working
  Top8rowidorder <- c("318", "98", "137", "142", "102", "432", "478", "112")
  
#Plot the antigen dilution curves for the best 8 conditions, 4 antigens (and GST?!?! or ALL antigens?!?!?!)
  
  #prepare transposed hello data frame
  hello2 <- tibble::column_to_rownames(hello, var = "rowid")
  helloT <- t(hello2) 
  colnames(helloT) <- c(as.character(rownames(hello2)))
  
  #make a target data frame with the new names because the current names don't match
  targetnames <- make.names(target_meta.df$Name)
  targethello <- target_meta.df
  targethello$Name <- targetnames
  
  #merge target metadata with data for the 8 conditions
  #Note: The print buffer column is meaningless at this point, ignore it!! Print buffer is part of the condition number
  hellocon <- merge(targethello, helloT, all.x = FALSE, all.y = FALSE, by.x = "Name", by.y = "row.names", sort = FALSE)
  
  helloconmelt <- melt(hellocon, measure.vars = c(colnames(helloT)))
  
  #FYI - getting a warning message: attributes are not identical across measure variables; they will be dropped,
  #but the data looks fine, 
  
  #FOR EACH ANTIGEN SEPARATELY - plot the concentration vs. condition
  aglist <- c("AMA1", "MSP1.19", "Hyp2", "EPF1v2", "GEXP18", "GST")
  
  #need to set value to numeric and concentration to factor, for some reason they aren't
  for(i in 1:length(aglist)){
    
    antigen <- aglist[[i]]
    
    agnums <- grep(antigen, helloconmelt$Name)
    oneantigen <- helloconmelt[agnums,]
    
    #explicitly factor level to match above
    oneantigen$variable <- factor(oneantigen$variable, levels = unique(Top8rowidorder))
    
    png(filename = paste0("CP3.",antigen, ".Dilutions.tif"), width = 4, height = 3.3, units = "in", res = 1200)
    
    print(ggplot(oneantigen, aes(x = as.factor(Concentration), y = as.numeric(value), color = variable)) + geom_point(shape=18, size = 2) +
      geom_line(aes(group = variable)) + theme_bw() + labs(y = "Normalized Log2(Positive/Negative)", x= "Concentration (µg/mL)", title = antigen) + 
      theme(axis.text.x = element_text(color = "black"), panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank()) +
      scale_color_hue(name = "Condition"))
    
    graphics.off()
    
  }
  
  remove(i)
  

#for all 480 conditions, CP3 only, plot the antigen dilutions versus log2(pos/neg), all antigens on one plot
  
  #test plot without specifying dilution and antigen - looks good, go ahead and do it properly
  CP3allmelt <- melt(CP3all)
  
  ggplot(CP3allmelt, aes(x = variable, y = value)) + geom_boxplot() +
          theme_bw() + labs(y = "Normalized Log2(Positive/Negative)", x= "Concentration/Antigen (?g/mL)") + 
          theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank()) +
          ylim(-2,10) + theme(axis.text.x = element_text(angle = 90, size = 8.5, vjust = 0.5, color = "black"))
  
  #base data is CP3all - need to link targets with dilutions - prepare transposed DF
  CP3all2 <- tibble::column_to_rownames(CP3all, var = "rowid")
  CP3allT <- t(CP3all2) 
  colnames(CP3allT) <- c(as.character(rownames(CP3all2)))
  
  #updated target_metadata is targethello
  CP3alltarget <- merge(targethello, CP3allT, all.x = FALSE, all.y = FALSE, by.x = "Name", by.y = "row.names", sort = FALSE)
  
  #need to add antigen labels to the data frame, this is column number 486
  CP3alltarget$Antigen <- c("AMA1")
  #use aglist created above to label antigens
  for(i in 1:length(aglist)){
    antigen <- aglist[[i]]
    nums <- grep(antigen, CP3alltarget$Name)
    CP3alltarget[nums, 486] <- antigen
  }
  
  #melt the data - getting Warning message: attributes are not identical across measure variables; they will be dropped, but data looks fine :/
  CP3antmelt <- melt(CP3alltarget, id.vars = c("Name", "Category", "Concentration", "Antigen", "Print_Buffer", "Expression_Tag"))
  
  #plot - this is all already GST subtracted data and is pos/neg ratio
  png(filename = paste0("CP3.ALL.Dilutions.tif"), width = 7, height = 3.5, units = "in", res = 1200)
  
  ggplot(CP3antmelt, aes(x = as.factor(Antigen), y = as.numeric(value), fill = as.factor(Concentration))) + geom_boxplot(outlier.size = 0.3) +
    theme_bw() + labs(y = "Normalized Log2(Positive/Negative)", x= "Antigen") + 
    theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank()) +
    theme(axis.text.x = element_text(color = "black")) + ylim(0,10) +
    scale_fill_hue(name = "Concentration (?g/mL)")
  
  graphics.off()
  
  #plot without GST
  CP3noGST <- filter(CP3antmelt, !(Antigen == "GST"))
  
  png(filename = paste0("CP3.noGST.Dilutions.tif"), width = 7, height = 3.5, units = "in", res = 1200)
  
  ggplot(CP3noGST, aes(x = as.factor(Antigen), y = as.numeric(value), fill = as.factor(Concentration))) + geom_boxplot(outlier.size = 0.3) +
    theme_bw() + labs(y = "Normalized Log2(Positive/Negative)", x= "Antigen") + 
    theme(panel.border = element_blank(), axis.line = element_line(), panel.grid = element_blank()) +
    theme(axis.text.x = element_text(color = "black")) + ylim(0,10) +
    scale_fill_hue(name = "Concentration (?g/mL)")
  
  graphics.off()
  
save.image(file = "PostOptimizationTop10v2.RData")  
  
  