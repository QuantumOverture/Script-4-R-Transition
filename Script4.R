require("TDA")

# Command Line Arguements
args <- commandArgs()
dataSet <- "horlings" #args[4]
homDim <- 2 #as.numeric(args[5])
partNum <- 101 #as.numeric(args[6])
epsIncr <- 0.01 #round(as.numeric(args[7]), 3)
action <- "sect" #args[8]

# Setup correct file paths
setwd(file.path(getwd(),".."))
dictFile <- paste(dataSet,"_",action,"_dict_",as.character(partNum),".txt",sep="")
dictPath <- file.path(getwd(),"Data",dataSet,action,dictFile)

dataFile <- paste(dataSet,"_data.txt",sep = "")
dataPath <- file.path(getwd(),"Data",dataSet,dataFile)

# Read data and dictionary and remove unneccessary rows
inputList <- read.table(dataPath, sep = '\t',header = FALSE) # ? What parameters do I give script 2 to 3 to get here?
inputList <- inputList[-2,]
dictList <- read.table(dictPath, sep = '\t',header = FALSE)
dictList <- dictList[-1,]

# Using dimension 2 as window
cloudDim = 2




MaxDIM = 1 # ?
MaxSCA = 5 # ? 


Diag <- ripsDiag(X = Circles, maxdimension = MaxDIM, maxscale =MaxSCA,
                 library = "Dionysus", printProgress = TRUE)

