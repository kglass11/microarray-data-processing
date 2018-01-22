#####################################
###READING IN YOUR MICROARRAY data###
#####################################

###Create a folder, into which you should copy this script, your .gpr files (named 'Slide 1.gpr', 'Slide 2.gpr' etc.), and your sample list csv file, called 'Sample list.csv".
# Your sample list file needs four columns. There names must be exactly as written here, though the order of the samples does not matter:
#1.slide_no
#2.sample_id
#3.block_rep_1
#4.block_rep_2

###Install any packages you may need for this script

require("gtools")

###Clear the environnment


### Set the working directory and path for the folder where .gpr files are

setwd("I:/Drakeley Group/Protein microarrays/Experiments/310817 Bijagos Islands/Screen 2")
path=("I:/Drakeley Group/Protein microarrays/Experiments/310817 Bijagos Islands/Screen 2")

###Identify all .gpr files in your set path

slide_ids <- list.files(path=path,pattern="*.gpr")

###Read in the contents of each of the .gpr files you just identified
#This simply combines the data in a list format
#If your gpr files differ in terms of the number of rows of information before the header, you need to edit this command
#i.e skip=32 means it will skip the fiorst 32 rows, and header=T means it will read the 33rd ow as a header, and the 34th as the first row of data
#Note: slide_no_temp should be the last slide added to the list
#Slide number is assigned based on the file name. As files are created without editing by the genepix software - this seems like the easiest way of identifying files.

slides_list <- list()
for(i in 1:length(slide_ids)) { 
  
  setwd(path)
  slides_list[[i]] <- read.table(slide_ids[i],skip=34,sep="\t",header=T)
  slide_no_temp <- substr(slide_ids[[i]],7,nchar(slide_ids[[i]])-4)
  slides_list[[i]][,"slide_no"] <- slide_no_temp
}
remove(i)

###Name your slides in this list, according to their file names

names(slides_list) <- slide_ids

###Bind all data from the slide data list (slides.list) into a single dataframe

slides_all.df <- c()

for(i in 1:length(slides_list)) { 
  
  slides_all.df <- rbind(slides_all.df,slides_list[[i]])
}
remove(i)

###Read in list of sample IDs

samples.df <- read.csv("Analysis sample list 2.csv", header=T)

###Create a vector listing all of your samples, in the order they appear in your samples_list file
#In our samples_list files, the sample ID is always in column two - which is why that column is picked out here

samples <- as.character(samples.df[1:nrow(samples.df), 2])

###Now, make a new sample variable that is unique (to avoid issues with multiple blanks, etc.)

samples_unique <- c(paste(as.character(samples.df[1:nrow(samples.df), 2]), rownames(samples.df), sep = "_"))
samples.df[,5] <- samples_unique
colnames(samples.df)[5] <- "sample_id_unique"

###Create vectors indicating the number of slides, blocks, and samples
#Slide and sample number are determined automatically from the data you input, whereas block number is manual in this instance

index_slide <- as.numeric(length(slides_list))
index_block <- 32
index_sample <- as.numeric(length(samples))

###Assign your sample_ids to each row of the combined slide data (slides_all.df)
#The order of data in your samples.df file is irrelevant, as long as each sample ID is correctly matched to its slide and block numbers

slides_all.df$sample_id_unique <- c()

for(i in 1:dim(slides_all.df)[1]){
  
  print(i)
  row_ite<-slides_all.df[i,]
  block_ite<-row_ite$Block
  slide_ite<-row_ite$slide_no
  sample_info_1<-samples.df[which(samples.df$slide_no==slide_ite),]
  match<-block_ite%in%sample_info_1[,"block_rep_1"]
  if(match==TRUE){
    value<-which(block_ite==sample_info_1[,"block_rep_1"])
  }else{
    value<-which(block_ite==sample_info_1[,"block_rep_2"])
  }
  sample_info_2<-sample_info_1[value,"sample_id_unique"]
  slides_all.df$Sample[i]<-as.character(sample_info_2)
}
remove(i, sample_info_1, sample_info_2, match, row_ite, block_ite, slide_ite)

###Write slides_all.df to a file to keep as a csv in your directory
write.csv(slides_all.df,file="slidesall_combinedGPR.csv", row.names=T)

### Make a spot annotations dataframe
#Column 42 is the sample id column
annotation_targets.df <- slides_all.df[which(slides_all.df[,42]==samples_unique[1]),1:4]
annotation_targets.df <- cbind(row.names(annotation_targets.df), annotation_targets.df)
colnames(annotation_targets.df)[1] <- "target_id_unique"
index_target <- as.numeric(length(annotation_targets.df$target_id_unique))
rownames(annotation_targets.df) <- c(paste(rownames(annotation_targets.df), annotation_targets.df$Name, annotation_targets.df$Block, sep = "_"))