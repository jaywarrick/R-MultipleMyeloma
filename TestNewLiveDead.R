library(data.table)
library(foreign)
base1 <- 500
base2 <- 1000
nm <- 2
nl <- 3
ny <- 3
n1 <- base1*ny*nl*nm
n2 <- base2*ny*nl*nm
n <- n1+n2
mu1 <- -1
mu2 <- 1
x1 <- rnorm(n1, mean=mu1, sd=0.6)
x2 <- rnorm(n2, mean=mu2, sd=0.6)
x <- sample(c(x1,x2), size=(n/2))
tempData <- data.frame(Id=rep(rep(1:(base1+base2),each=nl*ny),nm), Extra=5, Experiment='hello', Array.X=0, Array.Y=1:ny, Location=rep(1:nl, each=ny), Measurement=rep(c('R','G'), each=n/nm), Value=c(rep(500,n/2), sample(500*exp(x),n/2)))
tempData <- data.table(tempData)
tempData <- tempData[,lapply(.SD, as.factor)]
tempData$Value <- as.numeric(as.character(tempData$Value))
write.arff(tempData, '/Users/jaywarrick/Desktop/A Sandbox/madeUpData.arff')

hist(sample(exp(x),n/2))
hist(log(sample(exp(x),n/2)))
results <- analyzeLiveDead(logRatioThreshold=1.2, nClusters=2, locationDimension="", compiledTablePath='/Users/jaywarrick/Desktop/A Sandbox/madeUpData.arff', jexFolder='/Users/jaywarrick/Desktop/A Sandbox')

# logRatio = log(G/R) = x -> G/R = exp(x) -> G = R*exp(x) -> R = G/exp(x)

jexTempRFolder <- '/Users/jaywarrick/Documents/JEX/Example Database/temp/RScriptTempFolder'
library('foreign')
sourceGitHubFile('jaywarrick','R-General','master','.Rprofile')
sourceGitHubFile('jaywarrick','R-MultipleMyeloma','master','NewLiveDead.R')
results <- analyzeLiveDead(compiledTablePath='/Users/jaywarrick/Documents/JEX/Example Database/Example Dataset/Cell_x0_y0/File-FileFromLoren/x0_y0.arff', jexFolder=jexTempRFolder, logRatioThreshold=0.0, nClusters=3, locationDimension='Location')
