rm(list=ls())
library(foreign)
source('~/.Rprofile')
# Also export data with "Class" column
# Also export X Y data for contextual analysis.

# First read in the data
duh <- read.arff("/Volumes/JEX Cruncher/JEX Databases/Dominique/Mutant vs WT/Cell_x0_y0/File-Output ARFF Table/x0_y0.arff")
duh$Class <- "Mutant"
data <- duh
duh <- read.arff("/Volumes/JEX Cruncher/JEX Databases/Dominique/Mutant vs WT/Cell_x1_y0/File-Output ARFF Table/x1_y0.arff")
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

# Try a simple Naive Bayes Classifier
