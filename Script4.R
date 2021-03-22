require("TDA")
require("stringr")
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


Maximum_Distance_Between_Points <- function(cloud, cloudDim){
  MaxDist = -1
  for(Row in 1:nrow(cloud)){
    for(SecondRow in 1:nrow(cloud)){
      if(SecondRow == Row){
        next
      }
      TempDist = 0
      for(dim in 1:cloudDim){
        TempDist = TempDist + (cloud[Row,dim]-cloud[SecondRow,dim])**2
      }
      
      if(sqrt(TempDist) > MaxDist){
        MaxDist = sqrt(TempDist)
      }
    }
  }
  
  return(MaxDist/3)
}

TruncNum <- function(Float){
  # Truncate number, to make it match the given increment's decimal spots
  return(trunc(Float*(10^decimalSpotsOfIncr))/(10^decimalSpotsOfIncr))
}


############################Core Functions##########################




FormatBarcode <- function(Barcode){
  # Step 1 -> Truncate each Death and Birth Value after the given increment's decimal place.
  
  # Create two seperate data frames for each betti type
  TempBetti0 <- data.frame("dimension" = c(),"Birth" = c(),"Death" = c(),stringsAsFactors = FALSE)
  TempBetti1 <- data.frame("dimension" = c(),"Birth" = c(),"Death" = c(),stringsAsFactors = FALSE)
  
  # Iterate through original barcode and build up TempBetti0 and TempBetti1
  for(Line in 2:dim(Barcode)[1]){
    # Truncate numbers and add them to their respective increment groups
    CurBirth <- ceiling(TruncNum(BarCode[Line,"Birth"][["Birth"]])/epsIncr) * epsIncr
    CurDeath <- ceiling(TruncNum(BarCode[Line,"Death"][["Death"]])/epsIncr) * epsIncr
    
    # Skip if the birth and death value are the same
    if(CurBirth == CurDeath){
      next
    }
    
    # Create temp row that will be appended to the correct data frame
    TempRow <- data.frame("dimension" = c(BarCode[Line,"dimension"][["dimension"]]),"Birth" = c(CurBirth),"Death" =c(CurDeath))
    
    # Append to correct data frame
    if(BarCode[Line,"dimension"][["dimension"]] == 0){
      TempBetti0 <- rbind(TempBetti0, TempRow)
    }else{
      TempBetti1 <- rbind(TempBetti1, TempRow)
    }
    
  }
  
  # Step 2.a -> Sort and create new barcode 
  
  # Order Betti 0 data frame and add infinity/end indicator
  TempBetti0 <- TempBetti0[order(TempBetti0$Death),]
  TempBetti0 <- rbind(TempBetti0, data.frame("dimension" = c(0),"Birth" = c(0),"Death" =c("inf")))
  
  # Add and return Betti1 inclusion if Betti1 points appear duing analysis
  if(dim(TempBetti1)[1] != 0){
    TempBetti1 <- TempBetti1[order(TempBetti1$Birth),]
    return(rbind(TempBetti0,TempBetti1))
  }

  # otherwise return only the betti 0 data frame
  return(TempBetti0)
  
  

}

JagLineGenerateBetti0 <- function(FormattedBarcode, MaxDist){
  TempEps <- epsIncr
  Result <- c(length(FormattedBarcode$Death[FormattedBarcode$dimension == "0"]))
  while(TempEps < MaxDist){
    NumHitsForCurEps <- length(FormattedBarcode$Death[FormattedBarcode$Death == as.character(TempEps) & FormattedBarcode$dimension == "0"])
    Result <- c(Result,Result[length(Result)]-NumHitsForCurEps)
    TempEps <- TempEps + epsIncr
  }
  
  return(Result)
}

JagLineGenerateBetti1 <- function(FormattedBarcode, MaxDist){
  TempEps <- 0
  Result <- c()
  while(TempEps < MaxDist){
    NumBornForCurEps <- length(FormattedBarcode$Birth[FormattedBarcode$Birth == as.character(TempEps) & FormattedBarcode$dimension == "1"])
    NumDeadForCurEps <- length(FormattedBarcode$Death[FormattedBarcode$Death == as.character(TempEps) & FormattedBarcode$dimension == "1"])
    
    if(length(Result) == 0){
      Result <- c(0,(NumBornForCurEps+Result[length(Result)])-NumDeadForCurEps)
    }else{
      Result <- c(Result,(NumBornForCurEps+Result[length(Result)])-NumDeadForCurEps)
    }
    TempEps <- TempEps + epsIncr
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

# Result Path and directory
resultsPath <- file.path(getwd(),"Results",dataSet,action,paste(as.character(cloudDim),"D",sep = ""),"Homology",as.character(epsIncr))
dir.create(resultsPath,recursive = TRUE,showWarnings = FALSE)

for(i in 1:nrow(dictList)){
  chr = as.numeric(dictList[i,1])
  arm = dictList[i,2]
  
  beg = as.numeric(dictList[i,3])
  end = as.numeric(dictList[i,4])
  
  seg = as.numeric(dictList[i,6])

  # Barcode folder
  dir.create(file.path(getwd(),"Results",dataSet,action,paste(as.character(cloudDim),"D",sep = ""),"Homology",as.character(epsIncr),as.character(chr)),showWarnings = FALSE)
  
  # Jag file
  JagPathB0 <- file.path(getwd(),"Results",dataSet,action,paste(as.character(cloudDim),"D",sep = ""),"Homology",str_glue("B0_{cloudDim}D_{dataSet}_{action}_{chr}{arm}_seg{seg}.txt"))
  JagPathB1 <- file.path(getwd(),"Results",dataSet,action,paste(as.character(cloudDim),"D",sep = ""),"Homology",str_glue("B1_{cloudDim}D_{dataSet}_{action}_{chr}{arm}_seg{seg}.txt"))
  file.create(JagPathB0)
  file.create(JagPathB1)
  
  for(j in 1:nrow(inputList)){
    profile = build_individual_profile(inputList,j,beg,end)
    cloud = build_cloud(profile, cloudDim)
    
    MaxDist = Maximum_Distance_Between_Points(cloud, cloudDim)[["tempPoint"]]
    MaxDIM = 1 # 1 = Loops/holes
    MaxSCA = 100 #  
    
    
    Diag <- ripsDiag(X = cloud, maxdimension = MaxDIM, maxscale =MaxSCA,
                     library = "Dionysus", printProgress = TRUE)
    
    
    BarCode <- Diag[[1]]
    
    # Save Barcode
    file.create(file.path(getwd(),"Results",dataSet,action,paste(as.character(cloudDim),"D",sep = ""),"Homology",as.character(epsIncr),as.character(chr),str_glue("Inter_{cloudDim}D_hom{homDim}_{dataSet}_{chr}{arm}_pat{j}_seg{seg}.txt")),showWarnings = FALSE)
    FormattedBarcode <- FormatBarcode(BarCode)
    write.table(FormattedBarcode,file.path(getwd(),"Results",dataSet,action,paste(as.character(cloudDim),"D",sep = ""),"Homology",as.character(epsIncr),as.character(chr),str_glue("Inter_{cloudDim}D_hom{homDim}_{dataSet}_{chr}{arm}_pat{j}_seg{seg}.txt")),sep="\t",col.names = FALSE,row.names = FALSE,quote = FALSE)
    
    JagLineBetti0 <- JagLineGenerateBetti0(FormattedBarcode,MaxDist)
    write(paste(JagLineBetti0, collapse="\t"),JagPathB0,append = TRUE)
    
    JagLineBetti1 <- JagLineGenerateBetti1(FormattedBarcode,MaxDist)
    write(paste(JagLineBetti1, collapse="\t"),JagPathB1,append = TRUE)
  }
}
