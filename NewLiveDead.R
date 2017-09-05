# Load necessary libraries
library('foreign')
library('data.table')
library('EMCluster')

analyzeLiveDead <- function(compiledTablePath, jexFolder, logRatioThreshold=0, nClusters=2, seed=543, locationDimension='Location')
{
	# Define some helper functions
	plotHist <- function(data, LDClass, thresh, jexFolder, x, y, loc, prefix)
	{
		imagePath <- file.path(jexFolder,paste0(prefix, x, '_', y, '_', loc, '.tif'))
		tiff(filename=imagePath)
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
					hist(data[LDClass==i], breaks=myBreaks, main=paste0('X:', x, ' Y:', y, ' Loc:', loc), col=freshCol, add=add, freq=TRUE, ylim=c(0,myLim))
				}
				else
				{
					hist(data[LDClass==i], breaks=myBreaks, main=paste0('X:', x, ' Y:', y, ' Loc:', loc), col=freshCol, add=add, freq=TRUE, ylim=c(0,myLim))
				}
				add <- TRUE
			}
		}
		abline(v=thresh, lwd=2, col='blue')
		dev.off()
		return(imagePath)
	}

	plotHistAll <- function(data, LDClass, thresh, jexFolder, isLog, prefix)
	{
		imagePath <- file.path(jexFolder,paste0(prefix, '.tif'))
		tiff(filename=imagePath)
		if(isLog)
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
					hist(data[LDClass==i], breaks=myBreaks, main='log(Ratio) All', col=freshCol, add=add, freq=TRUE, ylim=c(0,myLim))
					add <- TRUE
				}
			}
		}
		else
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
					hist(data[LDClass==i], breaks=myBreaks, main='Ratio All', col=freshCol, add=add, freq=TRUE, ylim=c(0,myLim))
					add <- TRUE
				}
			}
		}
		abline(v=thresh, lwd=2, col='blue')
		dev.off()
		return(imagePath)
	}

	assignToClusters <- function(data, nClusters=2, rndSeed=1234)
	{
		set.seed(rndSeed)
		yo <- data[!is.na(data)]
		x <- data.frame(x=yo)


		emobj <- simple.init(x, nclass = nClusters)
		control <- .EMControl(alpha = 0.99, short.iter = 200, short.eps = 1e-2,
						  fixed.iter = 1, n.candidate = 3,
						  EM.iter = 100, EM.eps = 1e-3, exhaust.iter = 5)
		ret <- emcluster(x, emobj, assign.class = TRUE, EMC=control)

		temp <- data.frame(x=x$x, LDClass=as.numeric(as.character(ret$class)))
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

		temp$LDClass3 <- max(temp$LDClass)[1]
		temp$LDClass2 <- temp$LDClass3

		thresh <- list()
		for(i in (nrow(tempMu2)-1):1)
		{
			tempThresh <- min(temp[temp$LDClass == i & temp$x > tempMu2$mu[i],'x'])
			if(!length(tempThresh)==0 && !is.infinite(tempThresh[1]))
			{
				thresh[[i]] <- tempThresh[1]
				temp$LDClass2[temp$x < tempThresh] <- i-1
			}
		}

		return(list(data=temp, mu=tempMu2$mu, thresh=thresh, emclusterObj=ret))
	}

	# Read in the file
	duh <- read.arff(compiledTablePath)
	duh <- reorganizeFeatureTable(duh, specialNames=c(), convertToNumeric=FALSE)
	duh <- data.frame(duh)

	#Calculate the live/dead signal ratio
	duh$ratio <- (duh$G+101) / (duh$R+101) # this is to avoid division by 0 and assumes the background correction step adds back 100 intensity units to images.
	duh$logRatio <- log(duh$ratio)

	# Assign cells as either live (1) or dead (0)
	duh$LD <- 1
	duh$LD[duh$logRatio < logRatioThreshold] <- 0

	# Attempt EM Clustering to determine live / dead
	results <- assignToClusters(duh$logRatio, nClusters=nClusters, rndSeed=seed)
	duh$LDClass <- results$data$LDClass
	duh$LDClass2 <- results$data$LDClass2

	# Plot the data by array location and
	duh <- data.table(duh)
	if(locationDimension != '')
	{
		temp1 <- duh[ , plotHist(data=ratio, LDClass=LDClass2, jexFolder=jexFolder, thresh=exp(logRatioThreshold), x=Array.X[1], y=Array.Y[1], loc=.BY[[4]], prefix='hist1_'), by=c('Array.X', 'Array.Y', 'Experiment', locationDimension)]
		temp1b <- duh[ , plotHist(data=logRatio, LDClass=LDClass2, jexFolder=jexFolder, thresh=logRatioThreshold, x=Array.X[1], y=Array.Y[1], loc=.BY[[4]], prefix='hist2_'), by=c('Array.X', 'Array.Y', 'Experiment', locationDimension)]
	}
	else
	{
		temp1 <- duh[ , plotHist(data=ratio, LDClass=LDClass2, jexFolder=jexFolder, thresh=exp(logRatioThreshold), x=Array.X[1], y=Array.Y[1], loc='', prefix='hist1_'), by=c('Array.X', 'Array.Y', 'Experiment')]
		temp1b <- duh[ , plotHist(data=logRatio, LDClass=LDClass2, jexFolder=jexFolder, thresh=logRatioThreshold, x=Array.X[1], y=Array.Y[1], loc='', prefix='hist2_'), by=c('Array.X', 'Array.Y', 'Experiment')]
	}
	temp2 <- duh[ , plotHistAll(data=ratio, LDClass=LDClass2, jexFolder=jexFolder, thresh=exp(logRatioThreshold), isLog=FALSE, prefix='hist1_All'), ]
	temp2b <- duh[ , plotHistAll(data=logRatio, LDClass=LDClass2, jexFolder=jexFolder, thresh=logRatioThreshold, isLog=TRUE, prefix='hist2_All'), ]

	indRatioHistograms <- as.vector(temp1$V1)
	indLogRatioHistograms <- as.vector(temp1b$V1)
	overallRatioHistogram <- as.vector(temp2)
	overallLogRatioHistogram <- as.vector(temp2b)

	summaryTable <- data.frame()
	if(locationDimension != '')
	{
		summaryTable <- duh[ , list(LDRatio_Threshold=sum(LD)/length(LD), LDRatio_Cluster=length(LDClass2[LDClass2 > 0])/length(LDClass2)), by=c('Array.X','Array.Y','Experiment',locationDimension)]
	}
	else
	{
		summaryTable <- duh[ , list(LDRatio_Threshold=sum(LD)/length(LD), LDRatio_Cluster=length(LDClass2[LDClass2 > 0])/length(LDClass2)), by=c('Array.X','Array.Y','Experiment')]
	}

	path1 <- file.path(jexFolder,'SummaryTable.csv')
	write.csv(x=summaryTable,file=path1)
	path2 <- file.path(jexFolder,'SingleCellTable.csv')
	write.csv(x=duh,file=path2)
	summaryTable <- c(path1)
	singleCellTable <- c(path2)
	return(list(indRatioHistograms=indRatioHistograms, indLogRatioHistograms=indLogRatioHistograms, overallRatioHistogram=overallRatioHistogram, overallLogRatioHistogram=overallLogRatioHistogram, summaryTable=summaryTable, singleCellTable=singleCellTable, updatedTable=duh, clusterResults=results))
}

analyzeRatio <- function(compiledTablePath, outputFolder, logRatioThreshold=0, nClusters=2, seed=543, locationDimension='Location', Acol, Bcol, offset=NULL)
{

	# Define some helper functions
	plotHist <- function(data, cluster, thresh, outputFolder, x, y, loc, prefix)
	{
		imagePath <- file.path(outputFolder,paste0(prefix, x, '_', y, '_', loc, '.tif'))
		tiff(filename=imagePath)
		add <- FALSE
		col <- list(red=0, green=0, blue=0, alpha=0.5)
		tempHist <- hist(data, breaks=20, plot=FALSE)
		myBreaks <- tempHist$breaks
		myLim <- max(tempHist$counts)
		for(i in unique(cluster))
		{
			if(length(data[cluster==i])>0)
			{
				tempCol <- col
				tempCol[(i%%3)+1] <- 1
				freshCol <- do.call(rgb, tempCol)
				if(add)
				{
					rData <- range(data)
					hist(data[cluster==i], breaks=myBreaks, main=paste0('X:', x, ' Y:', y, ' Loc:', loc), col=freshCol, add=add, freq=TRUE, ylim=c(0,myLim))
				}
				else
				{
					hist(data[cluster==i], breaks=myBreaks, main=paste0('X:', x, ' Y:', y, ' Loc:', loc), col=freshCol, add=add, freq=TRUE, ylim=c(0,myLim))
				}
				add <- TRUE
			}
		}
		abline(v=thresh, lwd=2, col='blue')
		dev.off()
		return(imagePath)
	}

	plotHistAll <- function(data, cluster, thresh, outputFolder, isLog, prefix)
	{
		imagePath <- file.path(outputFolder,paste0(prefix, '.tif'))
		tiff(filename=imagePath)
		if(isLog)
		{
			add <- FALSE
			col <- list(red=0, green=0, blue=0, alpha=0.5)
			tempHist <- hist(data, breaks=40, plot=FALSE)
			myBreaks <- tempHist$breaks
			myLim <- max(tempHist$counts)
			for(i in unique(cluster))
			{
				if(length(data[cluster==i])>0)
				{
					tempCol <- col
					tempCol[(i%%3)+1] <- 1
					freshCol <- do.call(rgb, tempCol)
					hist(data[cluster==i], breaks=myBreaks, main='log(Ratio) All', col=freshCol, add=add, freq=TRUE, ylim=c(0,myLim))
					add <- TRUE
				}
			}
		}
		else
		{
			add <- FALSE
			col <- list(red=0, green=0, blue=0, alpha=0.5)
			tempHist <- hist(data, breaks=40, plot=FALSE)
			myBreaks <- tempHist$breaks
			myLim <- max(tempHist$counts)
			for(i in unique(cluster))
			{
				if(length(data[cluster==i])>0)
				{
					tempCol <- col
					tempCol[(i%%3)+1] <- 1
					freshCol <- do.call(rgb, tempCol)
					hist(data[cluster==i], breaks=myBreaks, main='Ratio All', col=freshCol, add=add, freq=TRUE, ylim=c(0,myLim))
					add <- TRUE
				}
			}
		}
		abline(v=thresh, lwd=2, col='blue')
		dev.off()
		return(imagePath)
	}

	assignToClusters <- function(data, nClusters=2, rndSeed=1234)
	{
		library(EMCluster)
		set.seed(rndSeed)
		yo <- data[!is.na(data)]
		x <- data.frame(x=yo)

		# Get basic cluster results (results are potentially out of order)
		emobj <- simple.init(x, nclass = nClusters)
		control <- .EMControl(alpha = 0.99, short.iter = 200, short.eps = 1e-2,
						  fixed.iter = 1, n.candidate = 3,
						  EM.iter = 100, EM.eps = 1e-3, exhaust.iter = 5)
		ret <- emcluster(x, emobj, assign.class = TRUE, EMC=control)

		# Create a data.frame to return
		temp <- data.frame(x=x$x, Cluster.Raw=as.numeric(as.character(ret$class)))
		tempMu <- data.frame(mu=as.vector(ret$Mu), Cluster.Raw=1:nrow(ret$Mu))

		# Order the mu table so we can go through sequentially and rename the clusters in ascending order
		tempMu <- tempMu[order(tempMu$mu),]

		# Create two copies so that you can use one as an original and another as an edited version
		# Originals will be without the '2' while news will be with the '2'
		temp2 <- temp
		tempMu2 <- tempMu
		for(i in 1:nrow(tempMu))
		{
			# Go through mu's in ascending order and assign the ascending order class
			temp2[temp$Cluster.Raw==tempMu$Cluster.Raw[i],'Cluster.Raw'] <- i
			# Also rename the clusters in the duplicate mu table
			tempMu2$Cluster.Raw[i] <- i
		}

		duh <- max(temp2$Cluster.Raw)[1]
		temp2$Cluster.Clean <- duh

		thresh <- list()
		# Go in reverse order from the max cluster number down to 1
		for(i in nrow(tempMu2):2)
		{
			# Find the value that discriminates between each pair of clusters
			tempThresh <- max(temp2[temp2$Cluster.Raw == i-1 & temp2$x < tempMu2$mu[i],'x'])
			if(!length(tempThresh)==0 && !is.infinite(tempThresh[1]))
			{
				# Then we found a threshold
				thresh[[i-1]] <- tempThresh[1]
				# Assign everything below that threshold to the next lowest cluster
				temp2$Cluster.Clean[temp2$x <= tempThresh[1]] <- i-1
			}
		}

		pi <- c()
		n <- nrow(temp2)
		for(i in 1:max(temp2$Cluster.Clean))
		{
			pi <- c(pi, sum(temp2$Cluster.Clean == i)/n)
		}

		return(list(data=temp2, mu=tempMu2$mu, thresh=thresh, emclusterObj=ret, pi=pi))
	}

	# Read in the file
	duh <- fread(compiledTablePath, header=T)
	# duh <- reorganize(duh, specialNames=c(), convertToNumeric=FALSE)
	duh <- data.frame(duh)

	# Check if Acol and Bcol exist in the data.table
	if(!(Acol %in% names(duh)))
	{
		stop(paste0('The A column, ', Acol, ', does not exist. Please provide a valid A column name. Names available are, ', paste0(names(duh), collapse=' '), '. Aborting.'))
	}
	if(!(Bcol %in% names(duh)))
	{
		stop(paste0('The B column, ', Bcol, ', does not exist. Please provide a valid B column name. Names available are, ', paste0(names(duh), collapse=' '), '. Aborting.'))
	}

	# Check to see if the offset if sufficient?
	Amin <- min(duh[[Acol]], na.rm=T)
	Amin.pos <- min(duh[[Acol]][duh[[Acol]] > 0])
	Bmin <- min(duh[[Bcol]], na.rm=T)
	Bmin.pos <- min(duh[[Bcol]][duh[[Bcol]] > 0])
	the.min <- min(c(Amin,Bmin), na.rm=T)

	offsetA <- 0
	if(Amin < 0)
	{
		# If there are negative values, then add an offset to everyone equal to -1.01 * min
		offsetA <- -1.1 * Amin
		warning(paste0('Adjusting A column by an offset to eliminate zeros and negative values. Offset = ', offsetA))
	}
	else if(Amin == 0)
	{
		# If there are zero values, add 1.1 * the min positive value
		offsetA <- 1.1 * Amin.pos
		warning(paste0('Adjusting A column by an offset to eliminate zeros and negative values. Offset = ', offsetA))
	}

	offsetB <- 0
	if(Amin < 0)
	{
		# If there are negative values, then add an offset to everyone equal to -1.01 * min
		offsetB <- -1.1 * Bmin
		warning(paste0('Adjusting B column by an offset to eliminate zeros and negative values. Offset = ', offsetB))
	}
	else if(Bmin == 0)
	{
		# If there are zero values, add 1.1 * the min positive value
		offsetB <- 1.1 * Bmin.pos
		warning(paste0('Adjusting B column by an offset to eliminate zeros and negative values. Offset = ', offsetB))
	}

	#Calculate the A/B signal ratio
	duh$ratio <- (duh[[Acol]] + offsetA) / (duh[[Bcol]] + offsetB)
	duh$logRatio <- log(duh$ratio)

	# Assign cells as either live (1) or dead (0) using threshold
	duh$thresh.class<- 1
	duh$thresh.class[duh$logRatio < logRatioThreshold] <- 0

	# Attempt EM Clustering to determine live / dead
	results <- assignToClusters(duh$logRatio, nClusters=nClusters, rndSeed=seed)
	duh$cluster.class <- results$data$Cluster.Clean

	# Plot the data by array location and
	duh <- data.table(duh)
	if(locationDimension != '')
	{
		temp1 <- duh[ , plotHist(data=ratio, cluster=cluster.class, outputFolder=outputFolder, thresh=exp(logRatioThreshold), x=Array.X[1], y=Array.Y[1], loc=.BY[[4]], prefix='hist1_'), by=c('Array.X', 'Array.Y', 'Experiment', locationDimension)]
		temp1b <- duh[ , plotHist(data=logRatio, cluster=cluster.class, outputFolder=outputFolder, thresh=logRatioThreshold, x=Array.X[1], y=Array.Y[1], loc=.BY[[4]], prefix='hist2_'), by=c('Array.X', 'Array.Y', 'Experiment', locationDimension)]
	}
	else
	{
		temp1 <- duh[ , plotHist(data=ratio, cluster=cluster.class, outputFolder=outputFolder, thresh=exp(logRatioThreshold), x=Array.X[1], y=Array.Y[1], loc='', prefix='hist1_'), by=c('Array.X', 'Array.Y', 'Experiment')]
		temp1b <- duh[ , plotHist(data=logRatio, cluster=cluster.class, outputFolder=outputFolder, thresh=logRatioThreshold, x=Array.X[1], y=Array.Y[1], loc='', prefix='hist2_'), by=c('Array.X', 'Array.Y', 'Experiment')]
	}
	temp2 <- duh[ , plotHistAll(data=ratio, cluster=cluster.class, outputFolder=outputFolder, thresh=exp(logRatioThreshold), isLog=FALSE, prefix='hist1_All'), ]
	temp2b <- duh[ , plotHistAll(data=logRatio, cluster=cluster.class, outputFolder=outputFolder, thresh=logRatioThreshold, isLog=TRUE, prefix='hist2_All'), ]

	indRatioHistograms <- as.vector(temp1$V1)
	indLogRatioHistograms <- as.vector(temp1b$V1)
	overallRatioHistogram <- as.vector(temp2)
	overallLogRatioHistogram <- as.vector(temp2b)

	summaryTable <- data.frame()

	# Write a function that returns percentages for each subpopulation

	if(locationDimension != '')
	{
		summaryTable <- duh[ , list(Ratio_Threshold=sum(thresh.class)/length(thresh.class), Ratio_Cluster=length(cluster.class[cluster.class > 0])/length(cluster.class)), by=c('Array.X','Array.Y','Experiment',locationDimension)]
	}
	else
	{
		summaryTable.cluster <- duh[, {
			tot = .N
			mu = mean(logRatio, na.rm=T)
			.SD[, .(mu, submu=mean(logRatio, na.rm=T), tot, subtot=.N), by='cluster.class']
			} , by=c('Array.X','Array.Y','Experiment')]
		setorder(summaryTable.cluster, Experiment, Array.X, Array.Y, cluster.class)
		summaryTable.cluster[, fraction := subtot/tot]

		summaryTable.manual <- duh[, {
			tot = .N
			mu = mean(logRatio, na.rm=T)
			.SD[, .(mu, submu=mean(logRatio, na.rm=T), tot, subtot=.N), by='thresh.class']
		} , by=c('Array.X','Array.Y','Experiment')]
		setorder(summaryTable.manual, Experiment, Array.X, Array.Y, thresh.class)
		summaryTable.manual[, fraction := subtot/tot]
	}

	path1 <- file.path(outputFolder,'SummaryTable.Cluster.csv')
	write.csv(x=summaryTable.cluster,file=path1)
	path2 <- file.path(outputFolder,'SummaryTable.Manual.csv')
	write.csv(x=summaryTable.manual,file=path2)
	path3 <- file.path(outputFolder,'SingleCellTable.csv')
	write.csv(x=duh,file=path3)
	summaryTable.cluster <- c(path1)
	summaryTable.manual <- c(path2)
	singleCellTable <- c(path2)
	return(list(indRatioHistograms=indRatioHistograms, indLogRatioHistograms=indLogRatioHistograms, overallRatioHistogram=overallRatioHistogram, overallLogRatioHistogram=overallLogRatioHistogram, summaryTable.cluster=summaryTable.cluster, summaryTable.manual=summaryTable.manual, singleCellTable=singleCellTable, updatedTable=duh, clusterResults=results))
}

# something potentially good to know for this file blah<-'b', blah2<-'a', temp[,mean(unlist(.SD[,'a',with=FALSE])),by=c(blah)]

# something potentially good to know for this file blah<-'b', blah2<-'a', temp[,mean(unlist(.SD[,'a',with=FALSE])),by=c(blah)]
