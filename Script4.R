require("TDA")
options(scipen=999)
setwd("C:/Users/ismai/Desktop/Script4Update/Research/TAaCGH")
# Command Line Arguements
args <- commandArgs()
dataSet <- "horlings" #args[4]
homDim <- 2 #as.numeric(args[5])
partNum <- 101 #as.numeric(args[6])
epsIncr <- 0.01 #round(as.numeric(args[7]), 3)
action <- "sect" #args[8]


decimalSpotsOfIncr <- nchar(strsplit(as.character(epsIncr),"[.]")[[1]][2])

#### TO DO LIST ####

# Add build_individual_profile for my simulated data

####################


#################IMPORTANT INTERNAL FUNCTIONS#######################

build_individual_profile <- function(CGHdata , patient , start , finish){
  tProfile <- CGHdata[patient,(start+1):(finish + 1)]
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
        if(i+j-1 <= ncol(CGHprofile)){
          tempPoint <- c( tempPoint,round(as.numeric(CGHprofile[ 1,(i+j-1)]),digits=6) )
        }else{
          tempPoint <- c( tempPoint,round(as.numeric(CGHprofile[ 1,( (i+j-1) %% ncol(CGHprofile) )]),digits=6) )
        }
      #tempPoint <- c( tempPoint,round(as.numeric(CGHprofile[ 1,( (i+j-1) %% ncol(CGHprofile) )]),digits=6) )
    }
  tCloud <- rbind(tCloud,tempPoint)
  tempPoint <- c()
  }
  
return(tCloud)
}

ClearAndSortBarCode <- function(BarCode,EpsLimit){
  Result <- c()#data.frame("Dimension"=c(1),"Birth"=c(),"Death"=c("-"))
  Dim0 <- c()
  Dim1 <- data.frame("Dimension"=c(1),"Birth"=c(1),"Death"=c(1))
  for(Row in 1:nrow(BarCode)){

    Dimension = as.numeric(BarCode[Row,"dimension"])
    Birth = round(as.numeric(BarCode[Row,"Birth"]),decimalSpotsOfIncr)
    Death = round(as.numeric(BarCode[Row,"Death"]),decimalSpotsOfIncr)
    if(Birth != Death){
      if(Dimension == 1){
        Dim1<- rbind(Dim1,c(Dimension,Birth,Death))

      }else{
        Dim0 <- c(Dim0,Death)
      }
      
    }
  
  }
  
  
  Dim1 <- Dim1[-1,]
  if(nrow(Dim1)>1){
  Dim1$"Birth" <- as.numeric(as.character(Dim1$"Birth"))
  Dim1 <- Dim1[order(Dim1$"Birth"),]
  }
  Dim0 <- sort(Dim0)
  
  # Betti 0 first
  for(num in Dim0){
    if(num != EpsLimit){
      Result <- rbind(Result,c(0,0,epsIncr * ceiling(num/epsIncr)))
    }else{
      Result <- rbind(Result,c(0,0,"inf"))
    }
  }
  # Betti 1 second
  if(nrow(Dim1)>1){
    for(row in 1:nrow(Dim1)){
      print("!")
      Result <- rbind(Result,c(1,Dim1[row,2],epsIncr * ceiling(Dim1[row,3]/epsIncr)))
    }
  }
  
  return(Result)
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
    MaxSCA = 100 #  
    
    
    Diag <- ripsDiag(X = cloud, maxdimension = MaxDIM, maxscale =MaxSCA,
                     library = "Dionysus", printProgress = TRUE)
    
    
    BarCode <- Diag[[1]]

    ClearedAndSortedBarCode <- ClearAndSortBarCode(BarCode,MaxSCA)
    # Start for betti 0
    print(ClearedAndSortedBarCode)
    Start <- 0
    # Make sure cloud graph is working properly
    
  }
}







