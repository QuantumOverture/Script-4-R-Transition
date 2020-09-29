require("TDA")

# Command Line Arguements
args <- commandArgs()
dataSet <- "horlings" #args[4]
homDim <- 2 #as.numeric(args[5])
partNum <- 101 #as.numeric(args[6])
epsIncr <- 0.01 #round(as.numeric(args[7]), 3)
action <- "sect" #args[8]


#### TO DO LIST ####

# Add build_individual_profile for my simulated data

####################


#################IMPORTANT INTERNAL FUNCTIONS#######################

build_individual_profile <- function(CGHdata , patient , start , finish){
  tProfile <- CGHdata[patient,start:(finish + 1)]
  return(tProfile)
}

build_cloud <- function(CGHprofile , dim){
  tempPoint <- c()
  tCloud <- c()
  # Go through the number of points for that chromosome arm
  for(i in 1:ncol(CGHprofile)){
    # Have as many coordinates as dimensions
    for(j in 1:dim){
        # This ensures proper wraparound.
        tempPoint <- c( tempPoint,round(as.numeric(CGHprofile[ 1,(i+j) %% ncol(CGHprofile) +1 ]),digits=6) )
    }
  tCloud <- rbind(tCloud,tempPoint)
  tempPoint <- c()
  }
  
return(tCloud)
}


####################################################################

# Setup correct file paths
setwd(file.path(getwd(),".."))
dictFile <- paste(dataSet,"_",action,"_dict_",as.character(partNum),".txt",sep="")
dictPath <- file.path(getwd(),"Data",dataSet,action,dictFile)

dataFile <- paste(dataSet,"_data.txt",sep = "")
dataPath <- file.path(getwd(),"Data",dataSet,dataFile)

# Read data and dictionary and remove unneccessary rows
inputList <- read.table(dataPath, sep = '\t',header = FALSE,stringsAsFactors=FALSE) 
inputList <- tail(inputList,-2)
dictList <- read.table(dictPath, sep = '\t',header = FALSE,stringsAsFactors=FALSE)
dictList <- dictList[-1,]

# Using dimension 2 as window
cloudDim <- 2

# Result Path
reasultsPath <- file.path(getwd(),"..","Results",dataSet,action,paste(as.character(cloudDim),"D",sep = ""),"Homology")


for(i in 1:nrow(dictList)){
  chr = as.numeric(dictList[i,1])
  arm = dictList[i,2]
  
  beg = as.numeric(dictList[i,3])
  end = as.numeric(dictList[i,4])
  
  seg = as.numeric(dictList[i,6])

  
  
  for(j in 1:nrow(inputList)){
    profile = build_individual_profile(inputList,j,beg,end)
    cloud = build_cloud(profile, cloudDim)
    MaxDIM = 1 # 1 = Loops/holes
    MaxSCA = 0.3 #  
    
    
    Diag <- ripsDiag(X = cloud, maxdimension = MaxDIM, maxscale =MaxSCA,
                     library = "Dionysus", printProgress = TRUE)
  }
}

plot(Diag[["diagram"]])


}




