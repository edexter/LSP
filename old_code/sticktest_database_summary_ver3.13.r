###################################################################
#THIS SCRIPT GENERATES SOME BASIC SUMMARIES FROM THE STICKTEST DATABASE
###################################################################
#Version 3.13 (10.06.2022)

#LOAD REQUIRED PACKAGES
library(reshape)
library(tidyr)

###################################################################
#USER DEFINED VARIABLES - EDIT AS NEEDED
###################################################################

#The file path on your local machine where the database is stored
databaseFilePath<-"C:/Users/ericd/Desktop/sticktestDatabase_1and2and3.csv"
df<-read.csv(databaseFilePath) #Load database

#The location of the directory where you want files to be saved
outputPath<- "C:/Users/ericd/Desktop/"

#The file name for the output files
BinaryOutputFileName <- "resistotypes_binary.csv"
QuantitativeOutputFileName <- "resistotypes_quantitative.csv"
FullOutputFileName <- "resistotypes_full.csv"

#The minimum frequency of stick required to call a clone "S"
threshold<-0.5

#The complexity of the output files. If FALSE, output will be 
#summarized by each host/parasite/stick site combination. If TRUE, 
#data from all stick sites will be collapsed into host/parasite 
#combinations
simpleOutput = TRUE

#The level of aggregation for summarized records. The default level is by host.
#This can be changed by altering the next line. A few examples are shown below.
#Attempting to aggregate by more than 3 variables is not supported at the moment, 
#but can be added if needed.
recordsBy<-c("host") #Default
#recordsBy<-c("host","researcher") #Summarize by host and researcher 
#recordsBy<-c("host","researcher","test_date") #Summarize by host,researcher, and date

#Create a "missing" ID for projects that have no ID listed
if( any(is.na(df[,c(1:8)])) ) warning('Missing values in project ID will be filled with "ID_missing" ')
df$projectID<-replace_na(df$projectID,"ID_missing")

#Optional filters. Change as needed. An example line is below.
#df<-df[which(df$projectID == "SWISSPONDPANEL"),] #Filter by project ID

###################################################################
#CORE FUNCTIONS - NO EDITING REQUIRED
###################################################################

##################################################################
#Count the total number of assays per host/parasite/site combination
##################################################################

#Get stick site names
#sites<-colnames(df[-1:-28]) #The old static way of subsetting stick variables
sites<-colnames(df[which(nchar(colnames(df))<3)]) #This is the new dynamic way

#Convert the data to long format
scores<-melt(data = df, id.vars = c(recordsBy,"parasite"), measure.vars = sites)

#Remove empty or unscored wells from the data
scores<-subset(scores[which(scores$variable != "O"),]) #Remove unscored results
#scores<-scores[1:1000,]

#Create an output grouping list based on user selected variables
attach(scores) #Note that this function won't work unless "scores" is attached
myList<-list()
for (i in 1:length(recordsBy)){
  myList[[i]]<-get(recordsBy[i])
  myList[[length(recordsBy)+1]]<-scores$parasite
  myList[[length(recordsBy)+2]]<-scores$variable
}

#Count the total number of assays per host/parasite/site
my.fun<-function(x){length(x[!is.na(x)])}
totals<-aggregate(scores$value, by = myList, FUN = my.fun )

#####################################
########################################
#Convert data back to wide format
#Note that this code can accomodate up to 5 user input aggregating terms
#totals<-totals[1:10000,] #debugging. delete this line

if (length(recordsBy) == 1) {
  totals<-cast(totals,Group.1~Group.2+Group.3,value = "x", sum)
}

if (length(recordsBy) == 2 ) {
  totals<-cast(totals,Group.1 + Group.2 ~ Group.3 + Group.4,value = "x", sum)
}

if (length(recordsBy) == 3 ) {
  totals<-cast(totals,Group.1 + Group.2 + Group.3 ~ Group.4 + Group.5,value = "x", sum)
}

#Now find the total number of positive sticks
#This follows the same comment from above about user input groupings
sticks<-aggregate(scores$value, by = myList, FUN = sum,na.rm = TRUE )

#sticks<-sticks[1:10000,] #debugging. keep this line commented

if (length(recordsBy) == 1) {
  sticks<-cast(sticks,Group.1~Group.2+Group.3,value = "x", sum)
}

if (length(recordsBy) == 2 ) {
  sticks<-cast(sticks,Group.1 + Group.2 ~ Group.3 + Group.4,value = "x", sum)
}

if (length(recordsBy) == 3 ) {
  sticks<-cast(sticks,Group.1 + Group.2 + Group.3 ~ Group.4 + Group.5,value = "x", sum)
}

#Replace NA's with zero
sticks[is.na(sticks)] <- 0

###################################################################
##Create simplified resistotypes
##################################################################
#This function collapses all stick sites into a single host/parasite
#combination, and renames columns as needed. This step will be 
#skipped if the variable "simpleOutput" is set to FALSE
#Note that this has gotten really long due to adding in multiple subsetting variables

if (length(recordsBy) == 1) {
  if(simpleOutput == TRUE){
    simpleSticks<-sticks[grepl('_X', colnames(sticks))]
    simpleTotals<-totals[grepl('_X', colnames(totals))]
    
    simpleSticks<- cbind(sticks[,1],simpleTotals - simpleSticks)
    simpleTotals<- cbind(totals[,1],simpleTotals)
    
    colnames(simpleSticks) <- sub("_X*", "", colnames(simpleSticks))
    colnames(simpleTotals) <- sub("_X*", "", colnames(simpleTotals))
    
    colnames(simpleSticks)[1]<-"host"
    colnames(simpleTotals)[1]<-"host"
    
    sticks <- simpleSticks
    totals <- simpleTotals
  }
}

if (length(recordsBy) == 2) {
  if(simpleOutput == TRUE){
    simpleSticks<-sticks[grepl('_X', colnames(sticks))]
    simpleTotals<-totals[grepl('_X', colnames(totals))]
    
    simpleSticks<- cbind(sticks[,1],sticks[,2],simpleTotals - simpleSticks)
    simpleTotals<- cbind(totals[,1],totals[,2],simpleTotals)
    
    colnames(simpleSticks) <- sub("_X*", "", colnames(simpleSticks))
    colnames(simpleTotals) <- sub("_X*", "", colnames(simpleTotals))
    
    colnames(simpleSticks)[1]<-recordsBy[1]
    colnames(simpleTotals)[1]<-recordsBy[1]
    
    colnames(simpleSticks)[2]<-recordsBy[2]
    colnames(simpleTotals)[2]<-recordsBy[2]
    
    sticks <- simpleSticks
    totals <- simpleTotals
  }
}

if (length(recordsBy) == 3) {
  if(simpleOutput == TRUE){
    simpleSticks<-sticks[grepl('_X', colnames(sticks))]
    simpleTotals<-totals[grepl('_X', colnames(totals))]
    
    simpleSticks<- cbind(sticks[,1],sticks[,2],sticks[,3],simpleTotals - simpleSticks)
    simpleTotals<- cbind(totals[,1],totals[,2],totals[,3],simpleTotals)
    
    colnames(simpleSticks) <- sub("_X*", "", colnames(simpleSticks))
    colnames(simpleTotals) <- sub("_X*", "", colnames(simpleTotals))
    
    colnames(simpleSticks)[1]<-recordsBy[1]
    colnames(simpleTotals)[1]<-recordsBy[1]
    
    colnames(simpleSticks)[2]<-recordsBy[2]
    colnames(simpleTotals)[2]<-recordsBy[2]
    
    colnames(simpleSticks)[3]<-recordsBy[3]
    colnames(simpleTotals)[3]<-recordsBy[3]
    
    sticks <- simpleSticks
    totals <- simpleTotals
  }
}
####################################################################
#CREATE QUANTITATIVE RESISTOTYPES
####################################################################
#Remove the negative stick columns from the dataframe
sticks<-sticks[!grepl('_X', colnames(sticks))]
totals<-totals[!grepl('_X', colnames(totals))]

quantScores<-sticks

if (length(recordsBy) == 1) {
  
  quantScores<-cbind(sticks[,1],(sticks[-c(1:length(recordsBy))] / totals[-c(1:length(recordsBy))]))
  colnames(quantScores)[1]<-recordsBy[1]
  numColumns<-ncol(quantScores)
  
}

if (length(recordsBy) == 2) {
  
  quantScores<-cbind(sticks[,1],sticks[,2],(sticks[-c(1:length(recordsBy))] / totals[-c(1:length(recordsBy))]))
  colnames(quantScores)[1]<-recordsBy[1]
  colnames(quantScores)[2]<-recordsBy[2]
  numColumns<-ncol(quantScores)
  
}

if (length(recordsBy) == 3) {
  
  quantScores<-cbind(sticks[,1],sticks[,2],sticks[,3],(sticks[-c(1:length(recordsBy))] / totals[-c(1:length(recordsBy))]))
  colnames(quantScores)[1]<-recordsBy[1]
  colnames(quantScores)[2]<-recordsBy[2]
  colnames(quantScores)[3]<-recordsBy[3]
  numColumns<-ncol(quantScores)
  
}

#Create Binary resistotypes
binaryResist<-quantScores

binaryResist[,(length(recordsBy)+1):numColumns]<-apply(binaryResist[,(length(recordsBy)+1):numColumns],2, function(x) {ifelse(x >= threshold, "S","R")})

if (length(recordsBy) == 1) {
  colnames(binaryResist)[1] <-recordsBy[1]
}

if (length(recordsBy) == 2) {
  colnames(binaryResist)[1] <-recordsBy[1]
  colnames(binaryResist)[2] <-recordsBy[2]
}

if (length(recordsBy) == 3) {
  colnames(binaryResist)[1] <-recordsBy[1]
  colnames(binaryResist)[2] <-recordsBy[2]
  colnames(binaryResist)[3] <-recordsBy[3]
  
}

#APPEND BINARY RESISTOTYPES WITH CONFIDENCE SCORES
####################################################################
confScores<-sticks

for (i in (length(recordsBy)+1):numColumns){
  confScores[i]<-paste(sticks[,i],totals[,i],sep="|")
}

allData <- sticks

for (i in (length(recordsBy)+1):numColumns){
  allData[i]<-paste(binaryResist[,i],confScores[,i],sep="     ")
}

detach(scores)
####################################################################
#Export files
write.csv(binaryResist, file = paste(outputPath,BinaryOutputFileName,sep=""),row.names = FALSE)
write.csv(allData, file = paste(outputPath,FullOutputFileName,sep=""),row.names = FALSE)
write.csv(quantScores, file = paste(outputPath,QuantitativeOutputFileName,sep=""),row.names = FALSE)

