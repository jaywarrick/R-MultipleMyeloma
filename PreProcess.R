rm(list=ls())
library(CRImage)
library(foreign)
source('~/.Rprofile')
# Also export data with "Class" column
# Also export X Y data for contextual analysis.

folder <- '/Users/jaywarrick/Documents/MMB/Grants/2015 - Wisconsin 2020/Miyamoto UW2020'
w <- 5
h <- 5

# First read in the data
duh <- read.arff("/Volumes/BeebeBig/Jay/JEX Databases/Dominique/Mutant vs WT/Cell_x0_y0/File-Output ARFF Table/x0_y0.arff")
duh$Class <- "Mutant"
data <- duh
duh <- read.arff("/Volumes/BeebeBig/Jay/JEX Databases/Dominique/Mutant vs WT/Cell_x1_y0/File-Output ARFF Table/x1_y0.arff")
duh$Class <- "WT"
data <- rbind(data, duh)

# Then reorganize the table
data <- reorganizeFeatureTable(data, convertToNumeric=FALSE)

# Remove/Keep any features that we don't/do want for now.
# Think about creating an X, Y, and "clumped" feature for contextual analysis.

# If keeping Haralick features, combine measures for each direction by averaging to make "rotationally invariant".
# Find all names with Horizontal in them
temp <- data
data <- temp

hNames <- names(data)[grep("Horizontal",names(data))]
vNames <- gsub("Horizontal", "Vertical", hNames)
dNames <- gsub("Horizontal", "Diagonal", hNames)
adNames <- gsub("Horizontal", "AntiDiagonal", hNames)
avgNames <- gsub("Horizontal", "Avg", hNames)

haralickNames <- data.frame(H=hNames, V=vNames, D=dNames, AD=adNames, avg=avgNames, stringsAsFactors=FALSE)
myfunc <- function(row, theNames)
{
     return(mean(row[,theNames$H] + row[,theNames$V] + row[,theNames$D] + row[,theNames$AD]))
}

for(i in 1:nrow(haralickNames))
{
     data[,haralickNames[i,5]] <- (data[,haralickNames[i,1]] + data[,haralickNames[i,2]] + data[,haralickNames[i,3]] + data[,haralickNames[i,4]])/4
     data <- data[,!(names(data) %in% as.character(haralickNames[i,1:4]))]
}

# Reposition "Class" column as last column
classes <- data$Class
data <- data[,!(names(data) %in% c("Class"))]
data$Class <- as.factor(classes)


write.arff(data, '/Users/jaywarrick/Desktop/TestOutput.arff')
data <- read.arff('/Users/jaywarrick/Desktop/TestOutput.arff')

MT <- subset(data, Class=='Mutant')
WT <- subset(data, Class=='WT')
hist(subset(data, Class=='WT' & firstorder.mean.650X705M < 550 & firstorder.mean.650X705M > 490)$firstorder.mean.650X705M, breaks=seq(490,550,1), add=FALSE, col='red', xlab='p65 Intensity', main='')
hist(subset(data, Class=='Mutant' & firstorder.mean.650X705M < 550 & firstorder.mean.650X705M > 490)$firstorder.mean.650X705M, breaks=seq(490,550,1), add=TRUE, col='blue')

pdf(file=file.path(folder,'Area.pdf'), w=w, h=h)
hist(WT$geometric.Area.NA, add=FALSE, breaks=seq(0,100000,300), col='red', xlab='Area', main='', xlim=c(0,5000), ann=F)
hist(MT$geometric.Area.NA, add=TRUE, breaks=seq(0,100000,300), col='blue')
dev.off()


pdf(file=file.path(folder,'Area.pdf'), w=w, h=h)
hist(WT$imagemoment.NormalizedCentralMoment12.485X525M, add=FALSE, breaks = 400, col='red', xlab='Area', main='', ann=F)
hist(MT$imagemoment.NormalizedCentralMoment12.485X525M, add=TRUE, breaks = 400, col='blue')
dev.off()

t.test(WT$geometric.Area.NA, MT$geometric.Area.NA)

x1 <- MT$imagemoment.CentralMoment11.650X705M/1e15 + rnorm(MT$imagemoment.CentralMoment11.650X705M, mean=0, sd=0.03)
y1 <- MT$imagemoment.HuMoment3.485X525M * rnorm(MT$imagemoment.HuMoment3.485X525M, mean=1, sd=0.12)
x2 <- WT$imagemoment.CentralMoment11.650X705M/1e15 + rnorm(WT$imagemoment.CentralMoment11.650X705M, mean=0, sd=0.03)
y2 <- WT$imagemoment.HuMoment3.485X525M * rnorm(WT$imagemoment.HuMoment3.485X525M, mean=1, sd=0.12)
pdf(file=file.path(folder,'CentralMoment11vsHuMomment3.pdf'), w=w, h=h)
plot(x1, y1, pch=21, cex=0.5, col=rgb(0,0,1,0.2), bg=rgb(0,0,1,0.2), xlim=c(-1.7,0.5), xlab='RELA (p65), Central Moment 11', ylab='CellTracker Green, Hu Moment 3', ann=F, log='y')
points(x2, y2, pch=21, cex=0.5, col=rgb(1,0,0,0.2), bg=rgb(1,0,0,0.2))
dev.off()

imagemoment.CentralMoment21.485X525M
x1 <- MT$imagemoment.HuMoment1.485X525M# * rnorm(MT$imagemoment.HuMoment1.650X705M, mean=1, sd=0.01)
y1 <- MT$imagemoment.HuMoment3.485X525M# * rnorm(MT$imagemoment.HuMoment3.485X525M, mean=1, sd=0.01)
x2 <- WT$imagemoment.HuMoment1.485X525M# * rnorm(WT$imagemoment.HuMoment1.650X705M, mean=1, sd=0.01)
y2 <- WT$imagemoment.HuMoment3.485X525M# * rnorm(WT$imagemoment.HuMoment3.485X525M, mean=1, sd=0.01)
pdf(file=file.path(folder,'HuMoment1vsHuMomment3.pdf'), w=w, h=h)
plot(x1, y1, pch=21, cex=0.5, col=rgb(0,0,1,0.2), bg=rgb(0,0,1,0.2), xlim=c(0,0.00005), ylim=c(1e-6, 1e-5), xlab='RELA (p65), Central Moment 11', ylab='CellTracker Green, Hu Moment 3', ann=F, log='y')
plot(x1, y1, pch=21, cex=0.5, col=rgb(0,0,1,0.2), bg=rgb(0,0,1,0.2),  xlab='RELA (p65), Central Moment 11', ylab='CellTracker Green, Hu Moment 3', ann=F, log='y')
points(x2, y2, pch=21, cex=0.5, col=rgb(1,0,0,0.2), bg=rgb(1,0,0,0.2))
dev.off()

x1 <- MT$imagemoment.HuMoment1.650X705M# * rnorm(MT$imagemoment.HuMoment1.650X705M, mean=1, sd=0.01)
y1 <- MT$imagemoment.HuMoment3.485X525M# * rnorm(MT$imagemoment.HuMoment3.485X525M, mean=1, sd=0.01)
x2 <- WT$imagemoment.HuMoment1.650X705M# * rnorm(WT$imagemoment.HuMoment1.650X705M, mean=1, sd=0.01)
y2 <- WT$imagemoment.HuMoment3.485X525M# * rnorm(WT$imagemoment.HuMoment3.485X525M, mean=1, sd=0.01)
pdf(file=file.path(folder,'HuMoment1vsHuMomment3.pdf'), w=w, h=h)
plot(x1, y1, pch=21, cex=0.5, col=rgb(0,0,1,0.2), bg=rgb(0,0,1,0.2), xlim=range(x1), xlab='RELA (p65), Central Moment 11', ylab='CellTracker Green, Hu Moment 3', ann=F, log='xy')
points(x2, y2, pch=21, cex=0.5, col=rgb(1,0,0,0.2), bg=rgb(1,0,0,0.2))
dev.off()

x1 <- MT$firstorder.mean.650X705M
y1 <- MT$firstorder.mean.485X525M
x2 <- WT$firstorder.mean.650X705M
y2 <- WT$firstorder.mean.485X525M
n1 <- length(x1)
n2 <- length(x2)
pdf(file=file.path(folder,'Intensity.pdf'), w=w, h=h)
plot(x2-500, y2-500, pch=21, cex=0.5, cex.axis=1, col=rgb(0,0,1,0.2), bg=rgb(0,0,1,0.2), xlim=c(1,75), ylim=c(1,1500), xlab='RELA (p65), Central Moment 11', ylab='CellTracker Green, Hu Moment 3', ann=F, log='xy')
points(x1-500, y1-500, pch=21, cex=0.5, cex.axis=1, col=rgb(1,0,0,0.4), bg=rgb(1,0,0,0.4))
dev.off()

# Blank plot
pdf(file=file.path(folder,'BlankPlot.pdf'), w=3.7, h=h)
plot(c(),c(), xlab='', ylab='', xlim=c(-0.5,4.5), ylim=c(0, 1.2), log='', xaxt='n')
axis(1,at=c(0,1,2,3,4), labels=c(0,100,1,10,100))
dev.off()

# Blank plot
pdf(file=file.path(folder,'BlankMLCurve.pdf'), w=3.7, h=h)
plot(c(),c(), xlab='', ylab='', xlim=c(0,50), ylim=c(0, 100), log='')
dev.off()

plot(x, y, pch=21, cex=0.5, col=rgb(0,0,1,0.2), bg=rgb(0,0,1,0.2), xlim=range(data$imagemoment.CentralMoment11.650X705M), xlab='RELA (p65), Central Moment 11', ylab='CellTracker Green, Hu Moment 3', xaxt='n')
points(WT$imagemoment.CentralMoment11.650X705M/1e15, WT$imagemoment.HuMoment3.485X525M, pch=21, cex=0.5, col=rgb(1,0,0,0.2), bg=rgb(1,0,0,0.2))
axis(1, at=c(-1e-5, 0, 1e5, 1e10, 1e15))
MTdia <- sqrt(MT$geometric.Area.NA/pi)*(6.45/40)*2
WTdia <- sqrt(WT$geometric.Area.NA/pi)*(6.45/40)*2
mad(MTdia)/median(MTdia)
mad(WTdia)/median(WTdia)
(median(MTdia)-median(WTdia))/(median(MTdia))

MTp65 <- MT$firstorder.mean.650X705M
WTp65 <- WT$firstorder.mean.650X705M
mad(MTp65)/median(MTp65)
mad(WTp65)/median(WTp65)
(median(WTp65)-median(MTp65))/(median(WTp65))

# Try a simple Naive Bayes Classifier
