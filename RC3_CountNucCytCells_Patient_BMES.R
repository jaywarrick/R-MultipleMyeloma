rm(list=ls())
library(foreign)
library(data.table)
source('/Users/jaywarrick/Public/DropBox/GitHub/R-General/.Rprofile')

# getTable <- function(ds, x, y, name)
# {
# 	# Read in the cell data
# 	jData <- getData('Y:/Jay/JEXDatabases/Randomized Co-culture', ds=ds, x=x, y=y, type='File', name=name)
# 	if(is.null(jData))
# 	{
# 		return(NULL)
# 	}
# 	temp <- data.table(read.arff(jData$fileList[1]))
# 	temp <- reorganize(temp)
# 	setnames(temp, c('390X440','485X525'), c('Nuc','Cyt')) # Rename channels
# 	temp$Time <- gsub('hr','',temp$Time) # Allow making Time numeric
# 	temp <- as.data.table(lapply(temp, function(temp){as.numeric(as.character(temp))})) # Make all columns numeric
# 	return(temp)
# }
#
# getTables <- function(DSs, Xs, Ys, name, ...)
# {
# 	dtList <- list()
# 	for(ds in DSs)
# 	{
# 		for(x in Xs)
# 		{
# 			for(y in Ys)
# 			{
# 				temp <- getTable(ds=ds, x=x, y=y, name=name)
# 				if(!is.null(temp))
# 				{
# 					temp$ds <- ds
# 					temp$x <- x
# 					temp$y <- y
# 					dtList[[paste(ds, x, y, sep='.')]] <- temp
# 				}
# 			}
# 		}
# 	}
# 	ret <- rbindlist(dtList, use.names=T, ...)
# }
#
# # Get the cell information
# x <- getTables(DSs=c('20160913','20160915'), Xs=0:3, Ys=0:1, name='Data Table')
# fwrite(x, file.path = 'Y:/Jay/JEXDatabases/R Analysis/RC3 20161004/x_BMES.csv')
x <- fread(input = '/Volumes/Seagate Backup Plus Drive/Materials for BMES/R Analysis/RC3 20161004/x_Patient_BMES.csv')

# # Get the microwell maps information
# maps <- getTables(DSs=c('20160913','20160915'), Xs=0:3, Ys=0:1, name='Microwell Map')
# setnames(maps, c('Point','Region'), c('Region_16','Region_0')) # Point represents the region number at time 16hr while Region represents the region number at time 0 hr.
# fwrite(maps, file.path = 'Y:/Jay/JEXDatabases/R Analysis/RC3 20161004/maps_BMES.csv')
maps <- fread(input = '/Volumes/Seagate Backup Plus Drive/Materials for BMES/R Analysis/RC3 20161004/maps_Patient_BMES.csv')

# Maps now contains a full list of all the relevant Loc-Region numbers at time 0 in "Region_0"
# So, remove items in the cell table that don't have corresponding File-Loc-Region complex Ids
# First make the complex region Ids (cRegIds)
x$cRegId <- paste(x$ds, x$x, x$y, x$Loc, x$Region, sep='.')
maps$cRegId <- paste(maps$ds, maps$x, maps$y, maps$Loc, maps$Region_0, sep='.')
maps$cRegId_16 <- paste(maps$ds, maps$x, maps$y, maps$Loc, maps$Region_16, sep='.')

# First, rename regions at
x[Time==16]$cRegId <- maps$cRegId[match(x[Time==16]$cRegId, maps$cRegId_16)]

# Now that it is renamed, we can select by maps$cRegId
x <- x[cRegId %in% maps$cRegId]

getLR <- function(a, b)
{
     a[a < 1] <- 1
     b[b < 1] <- 1
     c <- log(a/b)
     return(c)
}

x[,MMLR:=getLR(MM,Nuc)]
x[,StromaLR:=getLR(Stroma,Nuc)]

# Hopefully no row is NA. If there are, the offset likely needs to be increase for calculating the LR
lapply(x,function(x){which(is.na(x))})

# Anything with more Cyt than Nuc is considered a live MM cell.
x$MMCell <- 0
x$MMCell <- as.numeric(x$MMLR > 0)
x$Stroma <- 0
x$Stroma <- as.numeric(x$StromaLR > 0)
x2 <- x[,list(Tot=sum(MMCell)+sum(Stroma), LiveMM=sum(MMCell), LiveStroma=sum(Stroma)),by=c('Loc','Region','Time','ds','x','y')]



# Reorganize the table to show values by time. Fill with zeros any regions that had cells at time 0 but not at time 16
x2 <- reorganize(x2, measurementCols=c('Time'), valueCols=c('Tot','LiveMM','LiveStroma'), fill=0)

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
write.csv(propTestResults, file='/Volumes/Seagate Backup Plus Drive/Materials for BMES/R Analysis/RC3 20161004/propTestResults_Patient_BMES.csv')

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
     ret <- dt[Type %in% c(p1, p2), list(p1=p1, p2=p2, p1x=x[Type==p1], p1n=n[Type==p1], p2x=x[Type==p2], p2n=n[Type==p2], p1Mean=mean[Type==p1], p2Mean=mean[Type==p2], change=getSign(mean[Type==p1], mean[Type==p2]), p.value=prop.test(x=c(x[Type==p1], x[Type==p2]), n=c(n[Type==p1], n[Type==p2]))[['p.value']])]
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


# Special Test
p1 <- 'Tot.3_BCult.2_Tx.0_MCult.2'
p2 <- 'Tot.3_BCult.2_Tx.1_MCult.2'
temp <- x2[,list(x=sum(LiveMM_16), n=.N*.BY$LiveMM_0), by=c('Type', 'LiveMM_0')]
temp <- temp[order(Type)]
p1 <- temp[Type == 'Tot.3_BCult.2_Tx.1_MCult.2' & LiveMM_0 == 2]
p2 <- temp[Type == 'Tot.3_BCult.2_Tx.1_MCult.2' & LiveMM_0 == 1]
prop.test(x=c(p1$x, p2$x), n=c(p1$n, p2$n))

thingsToTest <- data.frame(                    p1='Tot.1_BCult.1_Tx.0_MCult.1', p2='Tot.1_BCult.1_Tx.1_MCult.1', note = 'TxVsCtrl_MonoFromMono')
thingsToTest <- rbind(thingsToTest, data.frame(p1='Tot.2_BCult.1_Tx.0_MCult.1', p2='Tot.2_BCult.1_Tx.1_MCult.1', note = 'TxVsCtrl_MonoFromMono'))
thingsToTest <- rbind(thingsToTest, data.frame(p1='Tot.3_BCult.1_Tx.0_MCult.1', p2='Tot.3_BCult.1_Tx.1_MCult.1', note = 'TxVsCtrl_MonoFromMono'))
thingsToTest <- rbind(thingsToTest, data.frame(p1='Tot.4_BCult.1_Tx.0_MCult.1', p2='Tot.4_BCult.1_Tx.1_MCult.1', note = 'TxVsCtrl_MonoFromMono'))
thingsToTest <- rbind(thingsToTest, data.frame(p1='Tot.5_BCult.1_Tx.0_MCult.1', p2='Tot.5_BCult.1_Tx.1_MCult.1', note = 'TxVsCtrl_MonoFromMono'))
thingsToTest$p1 <- as.character(thingsToTest$p1)
thingsToTest$p2 <- as.character(thingsToTest$p2)
TxVsCtrl_MonoFromMono <- summarizeAll(props, thingsToTest)
TxVsCtrl_MonoFromMono

thingsToTest <- data.frame(                    p1='Tot.1_BCult.2_Tx.0_MCult.1', p2='Tot.1_BCult.2_Tx.1_MCult.1', note = 'TxVsCtrl_MonoFromCo')
thingsToTest <- rbind(thingsToTest, data.frame(p1='Tot.2_BCult.2_Tx.0_MCult.1', p2='Tot.2_BCult.2_Tx.1_MCult.1', note = 'TxVsCtrl_MonoFromCo'))
thingsToTest <- rbind(thingsToTest, data.frame(p1='Tot.3_BCult.2_Tx.0_MCult.1', p2='Tot.3_BCult.2_Tx.1_MCult.1', note = 'TxVsCtrl_MonoFromCo'))
thingsToTest <- rbind(thingsToTest, data.frame(p1='Tot.4_BCult.2_Tx.0_MCult.1', p2='Tot.4_BCult.2_Tx.1_MCult.1', note = 'TxVsCtrl_MonoFromCo'))
thingsToTest <- rbind(thingsToTest, data.frame(p1='Tot.5_BCult.2_Tx.0_MCult.1', p2='Tot.5_BCult.2_Tx.1_MCult.1', note = 'TxVsCtrl_MonoFromCo'))
thingsToTest$p1 <- as.character(thingsToTest$p1)
thingsToTest$p2 <- as.character(thingsToTest$p2)
TxVsCtrl_MonoFromCo <- summarizeAll(props, thingsToTest)
TxVsCtrl_MonoFromCo

thingsToTest <- data.frame(                    p1='Tot.2_BCult.2_Tx.0_MCult.2', p2='Tot.2_BCult.2_Tx.1_MCult.2', note = 'TxVsCtrl_CoFromCo')
thingsToTest <- rbind(thingsToTest, data.frame(p1='Tot.3_BCult.2_Tx.0_MCult.2', p2='Tot.3_BCult.2_Tx.1_MCult.2', note = 'TxVsCtrl_CoFromCo'))
thingsToTest <- rbind(thingsToTest, data.frame(p1='Tot.4_BCult.2_Tx.0_MCult.2', p2='Tot.4_BCult.2_Tx.1_MCult.2', note = 'TxVsCtrl_CoFromCo'))
thingsToTest <- rbind(thingsToTest, data.frame(p1='Tot.5_BCult.2_Tx.0_MCult.2', p2='Tot.5_BCult.2_Tx.1_MCult.2', note = 'TxVsCtrl_CoFromCo'))
thingsToTest$p1 <- as.character(thingsToTest$p1)
thingsToTest$p2 <- as.character(thingsToTest$p2)
TxVsCtrl_CoFromCo <- summarizeAll(props, thingsToTest)
TxVsCtrl_CoFromCo

thingsToTest <- data.frame(                    p1='Tot.2_BCult.2_Tx.1_MCult.1', p2='Tot.2_BCult.2_Tx.1_MCult.2', note = 'MonoVsCo_Treated')
thingsToTest <- rbind(thingsToTest, data.frame(p1='Tot.3_BCult.2_Tx.1_MCult.1', p2='Tot.3_BCult.2_Tx.1_MCult.2', note = 'MonoVsCo_Treated'))
thingsToTest <- rbind(thingsToTest, data.frame(p1='Tot.4_BCult.2_Tx.1_MCult.1', p2='Tot.4_BCult.2_Tx.1_MCult.2', note = 'MonoVsCo_Treated'))
thingsToTest <- rbind(thingsToTest, data.frame(p1='Tot.5_BCult.2_Tx.1_MCult.1', p2='Tot.5_BCult.2_Tx.1_MCult.2', note = 'MonoVsCo_Treated'))
thingsToTest$p1 <- as.character(thingsToTest$p1)
thingsToTest$p2 <- as.character(thingsToTest$p2)
MonoVsCo_Treated <- summarizeAll(props, thingsToTest)
MonoVsCo_Treated

thingsToTest <- data.frame(                    p1='Tot.2_BCult.2_Tx.0_MCult.1', p2='Tot.2_BCult.2_Tx.0_MCult.2', note = 'MonoVsCo_Control')
thingsToTest <- rbind(thingsToTest, data.frame(p1='Tot.3_BCult.2_Tx.0_MCult.1', p2='Tot.3_BCult.2_Tx.0_MCult.2', note = 'MonoVsCo_Control'))
thingsToTest <- rbind(thingsToTest, data.frame(p1='Tot.4_BCult.2_Tx.0_MCult.1', p2='Tot.4_BCult.2_Tx.0_MCult.2', note = 'MonoVsCo_Control'))
thingsToTest <- rbind(thingsToTest, data.frame(p1='Tot.5_BCult.2_Tx.0_MCult.1', p2='Tot.5_BCult.2_Tx.0_MCult.2', note = 'MonoVsCo_Control'))
thingsToTest$p1 <- as.character(thingsToTest$p1)
thingsToTest$p2 <- as.character(thingsToTest$p2)
MonoVsCo_Control <- summarizeAll(props, thingsToTest)
MonoVsCo_Control

CompiledTable <- rbindlist(list(TxVsCtrl_MonoFromMono=TxVsCtrl_MonoFromMono, TxVsCtrl_MonoFromCo=TxVsCtrl_MonoFromCo, TxVsCtrl_CoFromCo=TxVsCtrl_CoFromCo, MonoVsCo_Treated=MonoVsCo_Treated, MonoVsCo_Control=MonoVsCo_Control))
write.csv(x=CompiledTable, file='/Volumes/Seagate Backup Plus Drive/Materials for BMES/R Analysis/RC3 20161004/CompiledTable_Patient_BMES.csv')