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
#stopped video at 20 minutes into it.
