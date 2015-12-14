# Load necessary libraries
library('foreign')
library('data.table')
library('EMCluster')


plotHist <- function(data, LDClass, xlab='log(Ratio)', thresh=0)
{
	add <- FALSE
	col <- list(red=0, green=0, blue=0, alpha=0.5)
	tempHist <- hist(data, breaks=20, plot=FALSE)
	myBreaks <- tempHist$breaks
	myLim <- max(tempHist$counts)
	for(i in unique(LDClass))
	{
		if(length(data[LDClass==i])>0)
		{
			tempCol <- col
			tempCol[(i%%3)+1] <- 1
			freshCol <- do.call(rgb, tempCol)
			if(add)
			{
				rData <- range(data)
				hist(data[LDClass==i], breaks=myBreaks, col=freshCol, add=add, freq=TRUE, ylim=c(0,myLim))
			}
			else
			{
				hist(data[LDClass==i], breaks=myBreaks, xlab=xlab, col=freshCol, add=add, freq=TRUE, ylim=c(0,myLim))
			}
			add <- TRUE
		}
	}
	abline(v=thresh, lwd=2, col='blue')
}

assignToClusters <- function(data, nClusters=2, rndSeed=1234)
{
	set.seed(rndSeed)
	yo <- data[!is.na(data)]
	if(length(yo) > 10000)
	{
		x <- data.frame(x=sample(yo, 10000))
	}
	else
	{
		x <- data.frame(x=yo)
	}
	
	
	emobj <- simple.init(x, nclass = nClusters)
	control <- .EMControl(alpha = 0.99, short.iter = 200, short.eps = 1e-2,
					  fixed.iter = 1, n.candidate = 3,
					  EM.iter = 100, EM.eps = 1e-3, exhaust.iter = 5)
	ret <- emcluster(x, emobj, assign.class = TRUE, EMC=control)
	
	temp <- data.frame(x=x$x, LDClass=ret$class)
	tempMu <- data.frame(mu=as.vector(ret$Mu), LDClass=1:nrow(ret$Mu))
	tempMu <- tempMu[order(tempMu$mu),]
	print(tempMu)
	temp2 <- temp
	tempMu2 <- tempMu
	for(i in 1:nrow(tempMu))
	{
		temp[temp2$LDClass==tempMu$LDClass[i],'LDClass'] <- i
		tempMu2$LDClass[i] <- i
	}
	
	temp$LDClass <- temp$LDClass-1
	
	thresh <- list()
	for(i in (nrow(tempMu)-1):1)
	{
		tempThresh <- min(temp[temp$LDClass == i,'x'])
		if(!length(tempThresh)==0)
		{
			thresh[[i]] <- tempThresh[1]
		}
	}
	
	return(list(data=temp, mu=tempMu$mu, thresh=thresh, emclusterObj=ret))
}

# # Read in the file
# duh <- read.arff(compiledTablePath)
# duh <- reorganizeFeatureTable(duh, specialNames=c(), convertToNumeric=FALSE)
# duh <- data.frame(duh)
# 
# #Calculate the live/dead signal ratio
# duh$ratio <- (duh$G+101) / (duh$R+101) # this is to avoid division by 0 and assumes the background correction step adds back 100 intensity units to images.
# duh$logRatio <- log(duh$ratio)

duh <- NULL
duh <- data.frame(logRatio=singleImg$Ratio)

# Assign cells as either live (1) or dead (0)
threshold <- 0
duh$LD <- 1
duh$LD[duh$logRatio < threshold] <- 0

# Attempt EM Clustering to determine live / dead
results <- assignToClusters(duh$logRatio, nClusters=3, rndSeed=1234)
duh$LDClass <- results$data$LDClass
plotHist(duh$logRatio, duh$LDClass)
ß