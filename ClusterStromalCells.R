#rm(list=ls())

source('D:/GitHub/R-General/.Rprofile')
source('D:/GitHub/R-Informatics-private/HuMoments.R')
source('D:/GitHub/R-Informatics-private/Zernike.R')
source('D:/GitHub/R-Cytoprofiling/PreProcessingHelpers.R')

# source('~/.Rprofile')
# source('~/Public/DropBox/GitHub/R-Informatics-private/HuMoments.R')
# source('/Users/jaywarrick/Public/DropBox/GitHub/R-Informatics-private/Zernike.R')
# source('~/Public/DropBox/GitHub/R-Cytoprofiling/PreProcessingHelpers.R')

#sourceGitHubFile(user='jaywarrick', repo='R-Cytoprofiling', branch='master', file='PreProcessingHelpers.R')
library(data.table)
library(foreign)
library(kernlab)

# A1,2,3 - WT
# A4,5,6 - NES
# A7,9,10 - WT SecOnly
# A8,11,12 - NES SecOnly
#
# B10,11,12 - WT
# B7,8,9 - NES
# B3,4,6 - WT SecOnly
# B1,2,5 - NES SecOnly

#dir1 <- '/Volumes/Seagate Backup Plus Drive/Data Tables/Dominique 1'
dir1 <- 'G:/Data Tables/MM Stromal Clustering/Dataset Name/File - Output CSV Table'
Pt600 <- c('x0_y0.csv', 'x1_y0.csv', 'x2_y0.csv', 'x3_y0.csv')
Pt602 <- c('x0_y2.csv', 'x1_y2.csv', 'x2_y2.csv', 'x3_y2.csv')

set.seed(1234)
n <- NULL
tableList <- list()
tableList <- append(tableList, getTableList(dir1, Pt600, 'Pt600', expt=1, sampleSize = n))
tableList <- append(tableList, getTableList(dir1, Pt602, 'Pt602', expt=1, sampleSize = n))
x1 <- rbindlist(tableList, use.names=TRUE)
x1[,Class:=NULL]


fixColNames(x1)

# Make things easier to peruse
setorder(x1, Pt, Expt, file, Id, Label, MaskChannel, Measurement, ImageChannel)
x1$Id <- as.character(x1$Id) # Avoid standardizing the Id

# save(x1, file='/Volumes/Seagate Backup Plus Drive/Data Tables/x1_sampled.Rdata')
# save(x1, file='G:/Data Tables/x1_sampled.Rdata')

# load(file='G:/Data Tables/x1_sampled.Rdata')
# set.seed(1234)
x1b <- fixLongTableStringsInCol(x1, 'Measurement')
#x1b <- removeMeasurementNamesContaining(x1b,'DNZ')
#x1b <- removeMeasurementNamesContaining(x1b,'NUCw')
x1b <- removeMeasurementNamesContaining(x1b,'ZernikeCircleX')
x1b <- removeMeasurementNamesContaining(x1b,'ZernikeCircleY')
x1b <- removeMeasurementNamesContaining(x1b,'DNZernikeInnerCircleX')
x1b <- removeMeasurementNamesContaining(x1b,'DNZernikeOuterCircleY')
x1b <- replaceSubStringInAllRowsOfCol(x1b, '395 X 455M', '390 X 440', 'ImageChannel')
x1b <- replaceSubStringInAllRowsOfCol(x1b, '485 X 525 M', '485 X 525', 'ImageChannel')
x1b <- replaceSubStringInAllRowsOfCol(x1b, '560 X 607 M', '560 X 607', 'ImageChannel')
x1b <- replaceSubStringInAllRowsOfCol(x1b, '650 X 705 M', '648 X 684', 'ImageChannel')
# x1b <- x1b[ImageChannel != '560 X 607']
#x1b <- x1b[MaskChannel == 'WholeCell']

# Get rid of cells that have some NA data
naData <- unique(x1b[!is.finite(Value)]$cId)
x1b <- x1b[!(cId %in% naData)]
goodData <- unique(x1b[Measurement == 'ZernikeMag11_NUCwFIXED']$cId)
x1b <- x1b[(cId %in% goodData)]

# Do some calculations
x2 <- intIntensityNormalizeCentralMoments(x1b)
x2 <- meanNormalizeZernikeMoments(x2)
x2 <- calculateHuMoments(x2)
x2 <- calculateZernikeDotProduct(x2)

# Tempororarily make the table wide to calculate averages of Haralick over the different directions
x2 <- getWideTable(x2)
x2 <- calculateRMSofHaralick(x2)
x2 <- removeExtraneousColumns(x2)

# Get our long table back and reorder
x3 <- getLongTableFromTemplate(x2, x1b)
setorder(x3, Pt, Expt, file, Id, Label, MaskChannel, Measurement, ImageChannel)


# Perform robust standardization (x-median)/mad (entertain idea of not applying to histogram bins)
x3 <- standardizeLongData(x3)

# Generate a table of differences between measures for each MaskChannel/ImageChannel/Measurement combination
x3 <- refactor(x3)
diffs <- calculateChannelDifferences(x3)

# Standardize the difference data
diffs <- standardizeLongData(diffs)

# Merge it with the original dataset, merging MaskChannel and ImageChannel into MeasurementName
x3 <- rbindlist(list(x3,diffs), use.names = TRUE)
x3$Measurement <- paste(x3$Measurement, x3$MaskChannel, x3$ImageChannel, sep='_')
x3[,MaskChannel:=NULL]
x3[,ImageChannel:=NULL]

# Get a wide table for machine learning and plotting
x3 <- getWideTable(x3)

# Fix naming issues introduced by merging MaskChannel and ImageChannel names with Measurement name
x3 <- fixColNames(x3)

# Perform final sorting of columns of data for easier perusing
x3 <- sortColsByName(x3)
# shinyData[,lapply(.SD, function(x){if(is.factor(x)){return(as.factor(x))}else{return(x)}})]
shinyData <- copy(x3)

# getNumericCols(shinyData)[as.logical(as.vector(shinyData[,lapply(.SD, function(x){length(which(!is.finite(x)))>0}),.SDcols=getNumericCols(shinyData)]))]

# Write the data for potential analysis outside R
#write.csv(shinyData, file='/Users/jaywarrick/Documents/MMB/Grants/2016 - RO1 Cytoprofiling/shinyData_Dom3.csv', row.names=FALSE)
write.csv(shinyData, file='G:/Data Tables/MMStroma.csv', row.names=FALSE)
shinyData <- fread('G:/Data Tables/MMStroma.csv')
# Look at the data
#browseShinyData()

# Need to remove names like Id, Label, ImRow, ImCol, Z,
dataToTest <- copy(shinyData)
dataToTest[, c('Id','Label','ImCol','ImRow','Z','Loc','file','Expt'):=NULL]

dataToTest <- removeColsContainingNames(dataToTest, names=c('THISwFIXED'))
dataToTest <- removeColsContainingNames(dataToTest, names=c('NUCwPADDEDNUC'))
dataToTest <- removeColsContainingNames(dataToTest, names=c('NUCwFIXED'))
dataToTest <- removeColsContainingNames(dataToTest, names=c('Histogram'))
dataToTest <- removeColsContainingNames(dataToTest, names=c('LBP','minus'))
dataToTest <- removeColsContainingNames(dataToTest, names=c('RHu'))
dataToTest <- removeColsContainingNames(dataToTest, names=c('DNZ'))
dataToTest <- removeColsContainingNames(dataToTest, names=c('ZernikeMag5'))
dataToTest <- removeColsContainingNames(dataToTest, names=c('ZernikeDot5'))
dataToTest <- removeColsContainingNames(dataToTest, names=c('ImageMoments.Central'))
removeColsWithInfiniteVals(dataToTest)

library(cluster)
library(Rtsne)
library(plotly)
set.seed(1235)
rtsne_out <- Rtsne(as.matrix(dataToTest[,-('cId'),with=F]))
pts <- data.table(rtsne_out$Y)
pts$cId <- dataToTest$cId
kc <- kmeans(rtsne_out$Y, 5)
pts$cluster <- kc$cluster

xlab <- list(title = 'X')
ylab <- list(title = 'Y')
Cluster <- pts$cluster
plot_ly(mode='markers', x=pts$V1, y=pts$V2, color=Cluster, text=pts$cId) %>%
layout(xaxis = xlab, yaxis = ylab, hovermode="closest", showlegend=FALSE)

#saveCropMontages(dataset='G:/JEX Databases/MM Stroma Clustering 1/Dataset Name', object='Image-Crop Montage', pts=pts)


saveCropMontages <- function(dataset='G:/JEX Databases/MM Stroma Clustering 1/Dataset Name', object='Image-Crop Montage', pts)
{
	clusterRange <- range(pts$cluster)
	clusters <- seq(clusterRange[1],clusterRange[2],1)
	pts[,copyThumbToFolder(dataset, object, cId, cluster, dest='G:/Data Tables/MM Stromal Clustering/Clusters'),by=cId]
}

copyThumbToFolder <- function(dataset, object, cId, clustId, dest)
{
	print(dataset)
	print(object)
	print(cId)
	print(clustId)
	print(dest)
	from <- getFilePath(dataset, object, cId)
	print(paste0('from=', from))
	to <- file.path(dest, clustId)
	print(to)
	dir.create(to)
	print(cId)
	fname <- gsub('.csv','',cId,fixed=T)
	fname <- gsub('.','_',fname,fixed=T)
	fnames <- unlist(strsplit(fname,'_',fixed=T))
	fname <- paste0(fnames[2],'_',fnames[3],'_Loc',fnames[4],'_Id',fnames[5])
	to <- file.path(to, paste0(fname,'.tif'))
	print(to)
	file.copy(from=from, to=to, overwrite=T)
}

#getFilePath(dataset='G:/JEX Databases/MM Stroma Clustering 1/Dataset Name', object='Image-Crop Montage', cId='1.x3_y2.csv.9.36')

getFilePath <- function(dataset, object, cId)
{
	print(cId)
	strings <- unlist(strsplit(cId, '.', fixed=T))
	xy <-  unlist(strsplit(strings[2], '_', fixed=T))
	x <-  unlist(strsplit(xy[1],'x',fixed=T))[2]
	y <-  unlist(strsplit(xy[2],'y',fixed=T))[2]
	loc <- paste0('Loc',strings[4])
	id <- paste0('Id',strings[5])

	xy <- paste0('x',x,'_y',y)
	objectPath <- file.path(dataset, paste0('Cell_',xy), object)
	fileList <- list.files(path=objectPath)
	for(f in fileList)
	{
		if(grepl(xy, f) && grepl(loc, f) && grepl(id, f))
		{
			print('Found a matching cell...')
			return(file.path(objectPath,f))
		}
	}
	return(NULL)
}
