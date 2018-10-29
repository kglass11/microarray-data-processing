#Hefin Rhys youtube tutorial

#PCA tutorial 

setwd("/Users/Katie/Desktop/R Programming - Coursera")
getwd()

library(ggplot2)

data(iris)
head(iris)
summary(iris)

#PCA can only accept continuous variables
#scaling is really important when variables are measured on different scales
myPr <- prcomp(iris[,-5], scale = TRUE)
myPr

plot(iris$Sepal.Length, iris$Sepal.Width)
plot(scale(iris$Sepal.Length), scale(iris$Sepal.Width))

#usually get one fewer components than variables entered
summary(myPr)

#percent of variance in the data explained by that component 
#1 and 2 explains 95% of the variability 
#1 explains 73% etc

#scree plot (sp?) shows variance(square of standard deviation) for each component
plot(myPr, type ="l")

#biplot, which includes eigenvectors for each variable
biplot(myPr, scale = 0)

#extract PCA output (item x in the list)
str(myPr)

#only get first 2 components
iris2 <- cbind(iris, myPr$x[,1:2])

#plot with stat ellipse which shows 95% confidence interval
ggplot(iris2, aes(x=PC1, y = PC2, color = Species, fill=Species)) + 
  stat_ellipse(geom = "polygon", col = "black", alpha = 0.5) +
  geom_point(shape = 21, color = "black")

#correlations between original variables and principle components 
cor(iris[,1:4], iris2[,6:7])

##### Clustering Analysis Tutorial (kmeans, hierarchical, density, models) 

#scale the data (get z scores, everything centered around mean of 0 and standard deviations away)
irisScaled <- scale(iris[,-5])

#######kmeans - k is a hyperparameter that is user specified
fitK <- kmeans(irisScaled, 3) #3 species
fitK

#look at within cluster sum of squares by cluster, the percent is an indicator
#of how much variability in the data is explained by assigning clusters
plot(iris, col=fitK$cluster)

#choosing K
k <- list()
for(i in 1:10){
  k[[i]] <- kmeans(irisScaled, i)
}

k

ss <- list()
for(i in 1:10) {
  ss[[i]] <- k[[i]]$betweenss/k[[i]]$totss
}

#usually where there is a shoulder in the data you should pick that number.
#for this data hard to decide between 2 and 3 clusters
plot(1:10, ss, type = "b", ylab = "betweenss/totss", xlab = "cluster (k)")

#look at how the clusters match clusters you can see in the data already
#could pick 2 or 3 clusters from this as well
for(i in 1:4){
  plot(iris, col = k[[i]]$cluster)
}

###### Hierarchical clustering

#starts by putting each point in it's own cluster, then merges clusters
#gives information about how clusters are similar to each other

#need to supply a distance matrix, which is the distance of every point to every other point

d <- dist(irisScaled)

#usually need to try different algorithms, ward.D2 pre-selected dunno why though
fitH <- hclust(d, "ward.D2")
plot(fitH)
rect.hclust(fitH, k = 3, border = "red")

clusters <- cutree(fitH, k = 3)
clusters
plot(iris, col = clusters)

####### model-based clustering 
#this package looks at many models and uses maximizing BIC to select model type and number of clusters
install.packages("mclust")
library(mclust)

fitM <- Mclust(irisScaled)
fitM
plot(fitM)
#select plots in the console - see BIC to see how the model was chosen

####### density based clustering 
install.packages("dbscan")
library(dbscan)

#"A point p is a core point if at least minPts points are within distance ε 
#(ε is the maximum radius of the neighborhood from p) of it (including p). 
#Those points are said to be directly reachable from p." - wikipedia, DBSCAN

#minpts often the number of dimensions (variables) in data + 1

#work out input for eps parameter with kNNdist - look for knee / elbow in data
#look for more directions on choosing these parameters in dbscan info
kNNdistplot(irisScaled, k=3)
abline(h= 0.7, col = "red", lty = 2)

fitD <- dbscan(irisScaled, eps = 0.7 , minPts = 5)
fitD
#noise points here are the points which do not fit in either cluster
plot(iris, col = fitD$cluster)

#as the number of dimensions increases, the number of observations needed
#to distinguish clusters grows exponentially. Therefore it is often recommended 
#to perform a dimension reduction technique such as principle component analysis
#first and then use the results of that to perform the clustering analysis

#recommneded to run several clustering algorithms and then decide which 
#best captures patterns in the data