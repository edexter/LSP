#Set working directory
setwd("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP")

#Load required packages
library(reshape2)
library(ggplot2)

#Load LSP data and format LSP data
LSP <- read.delim("data/LSP_depth_stats/LSP_depth_genotypes.txt", header=TRUE)
#LSP$LSP1 <- ifelse(grepl("1",LSP$state),1,0)
#LSP$LSP2 <- ifelse(grepl("2",LSP$state),1,0)
#LSP$LSP3 <- ifelse(grepl("3",LSP$state),1,0)

# Append the haplotype count per sample
LSP$LSP1count <- ifelse(LSP$state == "LSP1", 2,
                    ifelse(LSP$state == "1/2" | LSP$state == "1/3", 1,
                           0))

LSP$LSP2count <- ifelse(LSP$state == "LSP2", 2,
                        ifelse(LSP$state == "1/2" | LSP$state == "2/3", 1,
                               0))

LSP$LSP3count <- ifelse(LSP$state == "LSP3", 2,
                        ifelse(LSP$state == "1/3" | LSP$state == "2/3", 1,
                               0))

# Sum the haplotypes together to make sure there were no errors
LSP$ALLcount <- rowSums(LSP[,7:9])

#Simplify and reformat the data frame for plotting
df <- LSP
df$event <- substr(df$ID, 1, 1)

#Append sampling dates
df$month <- df$event
df$month <- ifelse(df$event == "D", "06/06/2019", df$month)
df$month <- ifelse(df$event == "E", "20/06/2019", df$month)
df$month <- ifelse(df$event == "F", "02/07/2019", df$month)
df$month <- ifelse(df$event == "G", "19/07/2019", df$month)
df$month <- ifelse(df$event == "I", "18/08/2019", df$month)
df$month <- ifelse(df$event == "J", "02/09/2019", df$month)
df$month <- as.Date(df$month, format = "%d/%m/%Y")

# Get total number of daphnia sequenced per date
table(df$month)

#Categorize samples as pre-epidemic or post-epidemic
df$period <- ifelse(df$event == "D" | df$event == "E", "before","after")

#Total number of haplotypes recovered in the before period
beforeSum <- sum(df$LSP1count[df$period == "before"],
                 df$LSP2count[df$period == "before"],
                 df$LSP3count[df$period == "before"])

#Total number of haplotypes recovered in the before period
afterSum <- sum(df$LSP1count[df$period == "after"],
                df$LSP2count[df$period == "after"],
                df$LSP3count[df$period == "after"])

#Verify that the number equals twice the number of samples (256 expected)
sum(beforeSum,afterSum) / 2

#Get the raw number and percent for each LSP before the epidemic
sum(df$LSP1count[df$period == "before"])
sum(df$LSP1count[df$period == "before"]) / sum(df$ALLcount[df$period == "before"])

sum(df$LSP2count[df$period == "before"])
sum(df$LSP2count[df$period == "before"]) / sum(df$ALLcount[df$period == "before"])

sum(df$LSP3count[df$period == "before"])
sum(df$LSP3count[df$period == "before"]) / sum(df$ALLcount[df$period == "before"])

# Make sure the amount add up
sum(df$LSP1count[df$period == "before"]) + 
  sum(df$LSP2count[df$period == "before"]) +
  sum(df$LSP3count[df$period == "before"])

#Get the raw number and percent for each LSP after the epidemic
sum(df$LSP1count[df$period == "after"])
sum(df$LSP1count[df$period == "after"]) / sum(df$ALLcount[df$period == "after"])

sum(df$LSP2count[df$period == "after"])
sum(df$LSP2count[df$period == "after"]) / sum(df$ALLcount[df$period == "after"])

sum(df$LSP3count[df$period == "after"])
sum(df$LSP3count[df$period == "after"]) / sum(df$ALLcount[df$period == "after"])

# Make sure the amount add up
sum(df$LSP1count[df$period == "after"]) + 
  sum(df$LSP2count[df$period == "after"]) +
  sum(df$LSP3count[df$period == "after"])

temp <- aggregate(df$LSP1count ~ df$period, FUN = sum)
temp$`df$LSP1count`
temp$`df$LSP1count` / beforeSum * 100

#Get the raw number and percent of LSP 2 pre-epidemic
temp <- aggregate(df$LSP2count ~ df$period, FUN = sum)
temp$`df$LSP2count`
temp$`df$LSP2count` / beforeSum * 100
  
#Load field data
field <- read.csv("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/data/eric_swisspond_data_2019_2020.csv")
field <- field[1:11,]
field$Density.sample.magna <- as.numeric(field$Density.sample.magna)
field$Infection.prevalence.... <- as.numeric(field$Infection.prevalence....)
field$Infection.prevalence....[c(1,2,11)] <- 0
field$date <- as.Date(field[,1],format = "%d.%m.%Y")

#Merge field and LSP data
df <- merge(field,df, by.x = "date", by.y = "month")

#Categorize samples as pre-epidemic or post-epidemic
df$period <- ifelse(df$event == "D" | df$event == "E", "before","after")

#LSP Barplot
df2 <- melt(df,id.vars = c("event","date","Infection.prevalence...."), measure.vars = c("LSP1","LSP2","LSP3"))
df2 <-df2[df2$value == 1,]
df2$infection <- as.numeric(df2$Infection.prevalence....)

df$LSP1count <- ifelse()

ggplot(df2, aes(x = date,group = variable, fill = variable)) + geom_bar(stat = "count", position = "dodge")+
  xlab("Sampling date") + ylab("Haplotype count")

ggplot(df2[], aes(x = date,group = variable, fill = variable)) + geom_bar(stat = "count", position = "fill")+
  xlab("Sampling date") + ylab("Haplotype count")

ggplot(df2, aes(x = event,y = infection)) + geom_bar(stat = "identity")

#Field data lineplot
ggplot(field, aes(x = date, y = Density.sample.magna))+geom_bar(stat = "identity")

df <- aggregate(df$value, by = list(df$month,df$variable), FUN = sum)
colnames(df) <- c("month","LSP","count")


#Standardize by 100
df2 <- df
totalD <- sum(df$count[df$month == "D"])
totalE <- sum(df$count[df$month == "E"])
totalF <- sum(df$count[df$month == "F"])
totalG <- sum(df$count[df$month == "G"])
totalI <- sum(df$count[df$month == "I"])
totalI <- sum(df$count[df$month == "J"])

cast(df)
ggplot(df2, aes(x = month, y = count, group = LSP, color = LSP)) + geom_line()


df2 <- df[df$value == 1,]
ggplot(df2[df2$month != "I",], aes(x = month,group = variable, fill = variable)) + geom_bar(stat = "count", position = "fill")


ggplot(df2, aes(x = month,group = LSP, color = LSP)) + geom_bar(stat = "count")
