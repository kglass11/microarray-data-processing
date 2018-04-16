#Optimization single antigen, single dilution plots 

antigen = "AMA1"
dilution = "100"
colnames(AHHHH.df[11])

 #all samples and all variables - Scatter plot faceted
 png(filename = paste0(antigen,"_", dilution, "_facet_grid.tif"), width = 8, height = 5, units = "in", res = 1200)
 par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))

ggplot(AHHHH.df, aes(x=slide_type, y=X13_1.PfAMA1.100ug.ml_1, color = sample, shape = blocking_buffer)) +
   geom_point(size = 1)  +
   theme_bw() +
   theme(text = element_text(size=9), axis.text.x = element_text(angle = 90, hjust=0.95, vjust = 0.1)) +
   labs(x = "Slide Type", y = "Normalized Log2(Positive/Negative)", title = paste(antigen, dilution, "µg/mL")) +
   facet_grid(print_buffer ~ block_dilution )

 graphics.off()

#violin plot of slide type and print buffer, everything else combined
 png(filename = paste0(antigen,"_", dilution, "_ST_PB.tif"), width = 8, height = 5, units = "in", res = 1200)
 par(mfrow=c(1,1), oma=c(3,1,1,1),mar=c(4.1,4.1,3.1,2.1))

ggplot(AHHHH.df, aes(x=slide_type, y=X13_1.PfAMA1.100ug.ml_1, color = print_buffer, fill = print_buffer)) +
   geom_violin(trim = FALSE) +
   theme_bw() +
   guides(color=FALSE) +
   theme(text = element_text(size=9), axis.text.x = element_text(angle = 45, hjust = 1), legend.position="top") +
   labs(x = "Slide Type", y = "Normalized Log2(Positive/Negative)", title = paste(antigen, dilution, "µg/mL"), fill = "Print Buffer")

graphics.off()

