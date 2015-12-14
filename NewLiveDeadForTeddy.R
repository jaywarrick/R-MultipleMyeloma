# Load necessary libraries
library('foreign')
library('data.table')
library('EMCluster')


plotHist <- function(data, LDClass, main=-1, xlab='log(Ratio)', thresh=0)
{
	add <- FALSE
	col <- list(red=0, green=0, blue=0, alpha=0.5)
	tempHist <- hist(data, breaks=40, plot=FALSE)
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
				hist(data[LDClass==i], breaks=myBreaks, xlab=xlab, main=main, col=freshCol, add=add, freq=TRUE, ylim=c(0,myLim))
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
	if(length(yo) > 50000)
	{
		x <- data.frame(x=sample(yo, 50000))
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
	
	return(list(data=temp, mu=tempMu2, thresh=thresh, emclusterObj=ret))
}

# Attempt EM Clustering to determine live / dead

##Change this section for every different file you use
bigTable <- read.csv('/Volumes/TeddyJEX/2015-12-08 MM Memory/MyExpt_Cells.csv')
bigTable$ratio <- bigTable$Intensity_MeanIntensity_Green/bigTable$Intensity_MeanIntensity_Red 


#Assign drug concentrations, where each 4 images goes 100, 10, 1, 0
bigTable$drugConcentration <- -1
bigTable$drugConcentration <- (bigTable$ImageNumber-1)%%4
bigTable$drugConcentration[bigTable$drugConcentration == 0] <- 100
bigTable$drugConcentration[bigTable$drugConcentration == 1] <- 10
bigTable$drugConcentration[bigTable$drugConcentration == 2] <- 1
bigTable$drugConcentration[bigTable$drugConcentration == 3] <- 0

#These are your days, they group by staining
bigTable$group <- -1
bigTable$group[bigTable$ImageNumber %in% c(1:16,49:56)] <- 1
bigTable$group[bigTable$ImageNumber %in% c(17:32)] <- 2
bigTable$group[bigTable$ImageNumber %in% c(33:48,57:80)] <- 3
bigTable$clusterNumber <- -1

#This part stays the same
for(i in unique(bigTable$group))
{
	badI = tryCatch({
	
	print(i)
	results <- assignToClusters(log(bigTable[bigTable$group == i,'ratio']), nClusters=2, rndSeed=1234)
	bigTable$clusterNumber[bigTable$group == i] <- results$data$LDClass
	print(results$mu)
	plotHist(results$data$x, results$data$LDClass, main=i, xlab = 'Log(Ratio)') },
	warning = function(w){
		print(paste0("This is crappy: ", i))
	}
	)
}


write.csv(x = bigTable , file = '/Volumes/TeddyJEX/2015-12-08 MM Memory/MyExpt_Cells.csv') 

### Trial codes

#days <- list(1=1:16, 2=17:)

