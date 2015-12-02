library(foreign)
library(data.table)
# duh <- read.arff('/Volumes/Data/JEX Databases/MM Samples/Pt564 40hrs/Cell_x0_y0/File-Counts (non-tumor, dead, compiled)/x0_y0.txt')
# duh$cell <- 'non-tumor'
# duh1 <- read.arff('/Volumes/Data/JEX Databases/MM Samples/Pt564 40hrs/Cell_x0_y0/File-Counts (non-tumor, live, compiled)/x0_y0.txt')
# duh1$cell <- 'non-tumor'
# duh <- rbind(duh, duh1)
# duh1 <- read.arff('/Volumes/Data/JEX Databases/MM Samples/Pt564 40hrs/Cell_x0_y0/File-Counts (tumor, dead, compiled)/x0_y0.txt')
# duh1$cell <- 'tumor'
# duh <- rbind(duh, duh1)
# duh1 <- read.arff('/Volumes/Data/JEX Databases/MM Samples/Pt564 40hrs/Cell_x0_y0/File-Counts (tumor, live, compiled)/x0_y0.txt')
# duh1$cell <- 'tumor'
# duh <- rbind(duh, duh1)

duh <- read.arff('/Volumes/Data/JEX Databases/MM Samples/Pt564 16hrs/Cell_x0_y0/File-Counts (non-tumor, dead, compiled)/x0_y0.txt')
duh$cell <- 'non-tumor'
duh1 <- read.arff('/Volumes/Data/JEX Databases/MM Samples/Pt564 16hrs/Cell_x0_y0/File-Counts (non-tumor, live, compiled)/x0_y0.txt')
duh1$cell <- 'non-tumor'
duh <- rbind(duh, duh1)
duh1 <- read.arff('/Volumes/Data/JEX Databases/MM Samples/Pt564 16hrs/Cell_x0_y0/File-Counts (tumor, dead, compiled)/x0_y0.txt')
duh1$cell <- 'tumor'
duh <- rbind(duh, duh1)
duh1 <- read.arff('/Volumes/Data/JEX Databases/MM Samples/Pt564 16hrs/Cell_x0_y0/File-Counts (tumor, live, compiled)/x0_y0.txt')
duh1$cell <- 'tumor'
duh <- rbind(duh, duh1)

duh$Channel <- sub("560 X 607","dead", duh$Channel)
duh$Channel <- sub("485 X 525","live", duh$Channel)

duh <- reorganizeFeatureTable(duh, specialNames=c(), convertToNumeric=FALSE, nameCol = 'Channel')

duh$fraction <- duh$live / (duh$live + duh$dead)

duh <- data.table(duh)
duh[,Mean:=mean(fraction),by=.(Tx,cell)]
duh[,p:=t.test(fraction[Tx=='Bort'],fraction[Tx=='Control'])$p.value, by=.(cell)]


jexDBFolder <- '/Volumes/Data/JEX Databases/MM Samples'
jexTempRFolder <- '/Volumes/Data/JEX Databases/MM Samples/temp/RScriptTempFolder'

library(foreign)
library(data.table)
source('~/.Rprofile')

data1 <- data.frame(name='FileName', Value='Pt564 16hrs/Cell_x0_y0/File-Counts (non-tumor, dead, compiled)/x0_y0.txt')
data2 <- data.frame(name='FileName', Value='Pt564 16hrs/Cell_x0_y0/File-Counts (non-tumor, live, compiled)/x0_y0.txt')
data3 <- data.frame(name='FileName', Value='Pt564 16hrs/Cell_x0_y0/File-Counts (tumor, dead, compiled)/x0_y0.txt')
data4 <- data.frame(name='FileName', Value='Pt564 16hrs/Cell_x0_y0/File-Counts (tumor, live, compiled)/x0_y0.txt')

duh <- read.arff(paste0(jexDBFolder,'/',data1$Value[1]))
duh$cell <- 'non-tumor'
duh1 <- read.arff(paste0(jexDBFolder,'/',data2$Value[1]))
duh1$cell <- 'non-tumor'
duh <- rbind(duh, duh1)
duh1 <- read.arff(paste0(jexDBFolder,'/',data3$Value[1]))
duh1$cell <- 'tumor'
duh <- rbind(duh, duh1)
duh1 <- read.arff(paste0(jexDBFolder,'/',data4$Value[1]))
duh1$cell <- 'tumor'
duh <- rbind(duh, duh1)

duh$Channel <- sub("560 X 607","dead", duh$Channel)
duh$Channel <- sub("485 X 525","live", duh$Channel)

duh <- reorganizeFeatureTable(duh, specialNames=c(), convertToNumeric=FALSE, nameCol = 'Channel')

duh$fraction <- duh$live / (duh$live + duh$dead)

duh <- data.table(duh)
duh[,Mean:=mean(fraction),by=.(Tx,cell)]
duh[,p:=t.test(fraction[Tx=='Bort'],fraction[Tx=='Control'])$p.value, by=.(cell)]

path <- paste0(jexTempRFolder,'/SummaryTable.csv')
write.csv(x=duh,file=path)
fileList1 <- c(path)
