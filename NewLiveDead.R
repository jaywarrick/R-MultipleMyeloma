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

          temp$LDClass2 <- max(temp$LDClass)[1]
          thresh <- list()
          for(i in (nrow(tempMu)-1):1)
          {
               tempThresh <- min(temp[temp$LDClass == i & temp$x > tempMu$mu[i-1],'x'])
               if(!length(tempThresh)==0)
               {
                    thresh[[i]] <- tempThresh[1]
                    temp$LDClass2[temp$x < tempThresh] <- i-1
               }
          }

          return(list(data=temp, mu=tempMu$mu, thresh=thresh, emclusterObj=ret))
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
          temp1b <- duh[ , plotHist(data=logRatio, LDClass=LDClass2, jexFolder=jexFolder, thresh=logRatioThreshold, LDClass=LDClass, x=Array.X[1], y=Array.Y[1], loc='', prefix='hist2_'), by=c('Array.X', 'Array.Y', 'Experiment')]
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

# something potentially good to know for this file blah<-'b', blah2<-'a', temp[,mean(unlist(.SD[,'a',with=FALSE])),by=c(blah)]