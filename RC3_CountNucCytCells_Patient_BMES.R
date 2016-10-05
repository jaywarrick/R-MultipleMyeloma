rm(list=ls())
library(foreign)
library(data.table)
source('D:/GitHub/R-General/.Rprofile')

getTable <- function(ds, x, y, name)
{
	# Read in the cell data
	jData <- getData('Y:/Jay/JEXDatabases/Randomized Co-culture', ds='20160913', x=x, y=y, type='File', name=name)
	if(is.null(jData))
	{
		return(NULL)
	}
	temp <- data.table(read.arff(jData$fileList[1]))
	temp <- reorganize(temp)
	setnames(temp, c('390X440','485X525'), c('Nuc','Cyt')) # Rename channels
	temp$Time <- gsub('hr','',temp$Time) # Allow making Time numeric
	temp <- as.data.table(lapply(temp, function(temp){as.numeric(as.character(temp))})) # Make all columns numeric
	return(temp)
}

getTables <- function(DSs, Xs, Ys, name, ...)
{
	dtList <- list()
	for(ds in DSs)
	{
		for(x in Xs)
		{
			for(y in Ys)
			{
				temp <- getTable(ds=ds, x=x, y=y, name=name)
				if(!is.null(temp))
				{
					temp$ds <- ds
					temp$x <- x
					temp$y <- y
					dtList[[paste(ds, x, y, sep='.')]] <- temp
				}
			}
		}
	}
	ret <- rbindlist(dtList, use.names=T, ...)
}

# Get the cell information
x <- getTables(DSs=c('20160913','20160916'), Xs=0:3, Ys=0:1, name='Data Table')
fwrite(x, file.path = 'Y:/Jay/JEXDatabases/R Analysis/RC3 20161004/x_Patient_BMES.csv')
x <- fread(input = 'Y:/Jay/JEXDatabases/R Analysis/RC3 20161004/x_Patient_BMES.csv')

# # Get the microwell maps information
maps <- getTables(DSs=c('20160913','20160916'), Xs=0:3, Ys=0:1, name='Microwell Map')
setnames(maps, c('Point','Region'), c('Region_16','Region_0')) # Point represents the region number at time 16hr while Region represents the region number at time 0 hr.
fwrite(maps, file.path = 'Y:/Jay/JEXDatabases/R Analysis/RC3 20161004/maps_Patient_BMES.csv')
maps <- fread(input = 'Y:/Jay/JEXDatabases/R Analysis/RC3 20161004/maps_Patient_BMES.csv')

# Maps now contains a full list of all the relevant Loc-Region numbers at time 0 in "Region_0"
# So, remove items in the cell table that don't have corresponding File-Loc-Region complex Ids
# First make the complex region Ids (cRegIds)
x$cRegId <- paste(x$ds, x$x, x$y, x$Loc, x$Region, sep='.')
maps$cRegId <- paste(x$ds, x$x, x$y, maps$Loc, maps$Region_0, sep='.')

# Rename Time 16 Regions to match those of Time 0 (had to do this because the image rotate too much... we needed
# different Regions for each timepoint).
x[Time == 16, Region:=maps$Region_0[match(Region, maps$Region_16)]]
x[Time == 16 | cRegId[Time == 0] %in% maps$cRegId] # First filter time 0 cRegIds according to Maps

##### THIS IS WHERE I NEED TO FOCUS #####
#  match(c(1,4,7),c(1,2,3,4,5,6,7)) >>> [1] 1 4 7

getLR <- function(a, b)
{
	a[a < 1] <- 1
	b[b < 1] <- 1
	c <- log(a/b)
	return(c)
}

x[,LR:=getLR(Cyt,Nuc)]

# Hopefully no row is NA. If there are, the offset likely needs to be increase for calculating the LR
lapply(x,function(x){which(is.na(x))})

# Anything with more Cyt than Nuc is considered a live MM cell.
x$MM <- 0
x$MM <- as.numeric(x$LR > 0)
x2 <- x[,list(Tot=length(MM), LiveMM=sum(MM)),by=c('Loc','Region','Time','ds','x','y')]

# Reorganize the table to show values by time. Fill with zeros any regions that had cells at time 0 but not at time 16
x2 <- reorganize(x2, measurementCols=c('Time'), valueCols=c('Tot','LiveMM'), fill=0)

# Only keep wells where there where there were cells to begin with.
x2 <- x2[Tot_0 > 0 & Tot_0 <= 5 & LiveMM_0 > 0 & (Tot_0 >= Tot_16)]


x2[,Culture:=(1 + as.numeric(LiveMM_0 < Tot_0)),]
x2[,Tx:=as.numeric(y==1),]
x2[,BECulture:=2-as.numeric(x==3),]
x2[,Type:=paste0('Tot.',Tot_0, '_BCult.', BECulture, '_Tx.', Tx, '_MCult.', Culture),]

n <- x2[, length(LiveMM_0), by='Type']
n <- n[order(n$Type)]
props <- x2[,list(x=sum(LiveMM_16), n=sum(LiveMM_0), mean=sum(LiveMM_16)/sum(LiveMM_0)), by='Type']
props <- props[order(Type)]
propTestResults <- as.data.frame(pairwise.prop.test(props$x, props$n, p.adjust.method = 'none')[[3]])
propTestResults <- cbind(data.frame(Type=props$Type[-1]), propTestResults)
names(propTestResults) <- c('Type', props$Type[-length(props$Type)])
fwrite(propTestResults, file.path='Y:/Jay/JEXDatabases/R Analysis/RC3 20161004/propTestResults_Patient_BMES.csv')

getSign <- function(m1, m2)
{
	temp <- m2-m1
	if(temp > 0)
	{
		return('+')
	}
	else if(temp == 0)
	{
		return('o')
	}
	else if(temp < 0)
	{
		return('-')
	}
}


summarize <- function(dt, p1, p2, note)
{
	if(sum(dt$Type %in% c(p1, p2)) != 2)
	{
		return(NULL)
	}
	ret <- dt[Type %in% c(p1, p2), list(p1=p1, p2=p2, p1x=x[Type==p1], p1n=n[Type==p1], p2x=x[Type==p2], p2n=n[Type==p2], p1Mean=mean[Type==p1], p2Mean=mean[Type==p2], change=getSign(dt$mean[dt$Type==p1], dt$mean[dt$Type==p2]), p.value=prop.test(x=c(x[Type==p1], x[Type==p2]), n=c(n[Type==p1], n[Type==p2]))[['p.value']])]
	ret$sig <- ''
	ret[p.value <= 0.05, sig:='*']
	ret[p.value <= 0.01, sig:='**']
	ret$note <- note
	return(ret)
}

summarizeAll <- function(dt, thingsToTest)
{
	dtList <- list()
	for(i in seq_along(thingsToTest$p1))
	{
		temp <- summarize(dt, thingsToTest$p1[i], thingsToTest$p2[i], thingsToTest$note[i])
		if(!is.null(temp))
		{
			dtList[[as.character(i)]] <- temp
		}
	}
	ret <- rbindlist(dtList, use.names=T)
}

thingsToTest <- data.frame(                    p1='Tot.1_BCult.1_Tx.0_MCult.1', p2='Tot.1_BCult.1_Tx.1_MCult.1', note = 'Drug effect in MonoBullsEye on 1 MM cells')
thingsToTest <- rbind(thingsToTest, data.frame(p1='Tot.2_BCult.1_Tx.0_MCult.1', p2='Tot.2_BCult.1_Tx.1_MCult.1', note = 'Drug effect in MonoBullsEye on 2 MM cells'))
thingsToTest <- rbind(thingsToTest, data.frame(p1='Tot.3_BCult.1_Tx.0_MCult.1', p2='Tot.3_BCult.1_Tx.1_MCult.1', note = 'Drug effect in MonoBullsEye on 3 MM cells'))
thingsToTest <- rbind(thingsToTest, data.frame(p1='Tot.4_BCult.1_Tx.0_MCult.1', p2='Tot.4_BCult.1_Tx.1_MCult.1', note = 'Drug effect in MonoBullsEye on 4 MM cells'))
thingsToTest <- rbind(thingsToTest, data.frame(p1='Tot.5_BCult.1_Tx.0_MCult.1', p2='Tot.5_BCult.1_Tx.1_MCult.1', note = 'Drug effect in MonoBullsEye on 5 MM cells'))
thingsToTest$p1 <- as.character(thingsToTest$p1)
thingsToTest$p2 <- as.character(thingsToTest$p2)
MonoDrugEffect <- summarizeAll(props, thingsToTest)
MonoDrugEffect

thingsToTest <- data.frame(                    p1='Tot.2_BCult.2_Tx.1_MCult.1', p2='Tot.2_BCult.2_Tx.1_MCult.2', note = 'Co-Culture effect on Drug Sensitivity of 2 cells in CoBullsEye.')
thingsToTest <- rbind(thingsToTest, data.frame(p1='Tot.3_BCult.2_Tx.1_MCult.1', p2='Tot.3_BCult.2_Tx.1_MCult.2', note = 'Co-Culture effect on Drug Sensitivity of 3 cells in CoBullsEye.'))
thingsToTest <- rbind(thingsToTest, data.frame(p1='Tot.4_BCult.2_Tx.1_MCult.1', p2='Tot.4_BCult.2_Tx.1_MCult.2', note = 'Co-Culture effect on Drug Sensitivity of 4 cells in CoBullsEye.'))
thingsToTest <- rbind(thingsToTest, data.frame(p1='Tot.5_BCult.2_Tx.1_MCult.1', p2='Tot.5_BCult.2_Tx.1_MCult.2', note = 'Co-Culture effect on Drug Sensitivity of 5 cells in CoBullsEye.'))
thingsToTest$p1 <- as.character(thingsToTest$p1)
thingsToTest$p2 <- as.character(thingsToTest$p2)
MonoVsCoDrugEffect <- summarizeAll(props, thingsToTest)
MonoVsCoDrugEffect

thingsToTest <- data.frame(                    p1='Tot.2_BCult.1_Tx.1_MCult.1', p2='Tot.2_BCult.2_Tx.1_MCult.2', note = 'Co-Culture effect on Drug Sensitivity of 2 cells in CoBullsEye.')
thingsToTest <- rbind(thingsToTest, data.frame(p1='Tot.3_BCult.1_Tx.1_MCult.1', p2='Tot.3_BCult.2_Tx.1_MCult.2', note = 'Co-Culture effect on Drug Sensitivity of 3 cells in CoBullsEye.'))
thingsToTest <- rbind(thingsToTest, data.frame(p1='Tot.4_BCult.1_Tx.1_MCult.1', p2='Tot.4_BCult.2_Tx.1_MCult.2', note = 'Co-Culture effect on Drug Sensitivity of 4 cells in CoBullsEye.'))
thingsToTest <- rbind(thingsToTest, data.frame(p1='Tot.5_BCult.1_Tx.1_MCult.1', p2='Tot.5_BCult.2_Tx.1_MCult.2', note = 'Co-Culture effect on Drug Sensitivity of 5 cells in CoBullsEye.'))
thingsToTest$p1 <- as.character(thingsToTest$p1)
thingsToTest$p2 <- as.character(thingsToTest$p2)
MonoVsCoDrugEffect2 <- summarizeAll(props, thingsToTest)
MonoVsCoDrugEffect2

thingsToTest <- data.frame(                    p1='Tot.2_BCult.2_Tx.0_MCult.1', p2='Tot.2_BCult.2_Tx.0_MCult.2', note = 'Co-Culture effect on Viability of 2 cells in CoBullsEye.')
thingsToTest <- rbind(thingsToTest, data.frame(p1='Tot.3_BCult.2_Tx.0_MCult.1', p2='Tot.3_BCult.2_Tx.0_MCult.2', note = 'Co-Culture effect on Viability of 3 cell in CoBullsEye.'))
thingsToTest <- rbind(thingsToTest, data.frame(p1='Tot.4_BCult.2_Tx.0_MCult.1', p2='Tot.4_BCult.2_Tx.0_MCult.2', note = 'Co-Culture effect on Viability of 4 cell in CoBullsEye.'))
thingsToTest <- rbind(thingsToTest, data.frame(p1='Tot.5_BCult.2_Tx.0_MCult.1', p2='Tot.5_BCult.2_Tx.0_MCult.2', note = 'Co-Culture effect on Viability of 5 cell in CoBullsEye.'))
thingsToTest$p1 <- as.character(thingsToTest$p1)
thingsToTest$p2 <- as.character(thingsToTest$p2)
MonoVsCoViability <- summarizeAll(props, thingsToTest)
MonoVsCoViability


# # Replace any new NAs with 0
# x4 <- reorganize(x3, idCols=c('Loc','Region'), measurementCols=c('Time'), valueCols=c('Tot','LiveMM'))
# x4$Culture <- 1 + as.numeric(x4$LiveMM_0 < x4$Tot_0) # 1 = Mono 2 = co-culture
# # Remove any times where the initial total number of cells was 0.
# 
# 
# # for(colname in names(x2))
# # {
# # 	print(colname)
# # 	hist(x2[[colname]], breaks=40, main=colname)
# # }
# # lapply(x2, min)
# # lapply(x2, hist, breaks=40, add=F)
# 
# # hist(x2[Time==0]$LR, breaks=40, xlim=c(-8,6))
# # hist(x2[Time==16]$LR, breaks=40, xlim=c(-8,6))
# # hist(log(x2[Time==0]$Cyt),breaks=40)
# # hist(log(x2[Time==16]$Cyt), breaks=40)
# 
# x <- reorganize(data=x, measurementCols=c('Measurement','Time'), valueCols='Value')
