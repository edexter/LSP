################################################################################
##MARKER 1 (Interlink)
#Load marker data
df <- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/markers/LSP_marker_1_genotypes_interlink.txt")

#Extract and format the genotypes
GT <- df[,grepl("GT", colnames(df))]
colnames(GT) <- gsub(".GT","",colnames(GT))
GT <- data.frame(t(GT))
GT$ID <-row.names(GT)
GT[,1] <- gsub("\\|","/",GT[,1])
colnames(GT)[1] <- "M1"

#Extract and format genotype qualities
GQ <- df[,grepl("GQ", colnames(df))]
colnames(GQ) <- gsub(".GQ","",colnames(GQ))
GQ <- data.frame(t(GQ))
GQ$ID <-row.names(GQ)
colnames(GQ)[1] <- "M1Q"

#Merge formatted genotype and quality data for marker
M1 <- merge(GT, GQ , by = "ID")

##MARKER 2
#Load marker data
df <- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/markers/LSP_marker_2_genotypes_interlink.txt")

#Extract and format the genotypes
GT <- df[,grepl("GT", colnames(df))]
colnames(GT) <- gsub(".GT","",colnames(GT))
GT <- data.frame(t(GT))
GT$ID <-row.names(GT)
GT[,1] <- gsub("\\|","/",GT[,1])
colnames(GT)[1] <- "M2"

#Extract and format genotype qualities
GQ <- df[,grepl("GQ", colnames(df))]
colnames(GQ) <- gsub(".GQ","",colnames(GQ))
GQ <- data.frame(t(GQ))
GQ$ID <-row.names(GQ)
colnames(GQ)[1] <- "M2Q"

#Merge formatted genotype and quality data for marker
M2 <- merge(GT, GQ , by = "ID")

##MARKER 3
#Load marker data
df <- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/markers/LSP_marker_3_genotypes_interlink.txt")

#Extract and format the genotypes
GT <- df[,grepl("GT", colnames(df))]
colnames(GT) <- gsub(".GT","",colnames(GT))
GT <- data.frame(t(GT))
GT$ID <-row.names(GT)
GT[,1] <- gsub("\\|","/",GT[,1])
colnames(GT)[1] <- "M3"

#Extract and format genotype qualities
GQ <- df[,grepl("GQ", colnames(df))]
colnames(GQ) <- gsub(".GQ","",colnames(GQ))
GQ <- data.frame(t(GQ))
GQ$ID <-row.names(GQ)
colnames(GQ)[1] <- "M3Q"

#Merge formatted genotype and quality data for marker
M3 <- merge(GT, GQ , by = "ID")

markers <- merge(M1, M2, by = "ID")
markers <- merge(markers, M3, by = "ID")

##LSP DEPTH
#Load LSP depth data
LSP <- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/data/LSP_depth_stats/LSP_depth_genotypes.txt")
markers <- merge(markers, LSP, by = "ID")

temp <- markers[markers$M1Q > 1 & markers$state != "Unknown",]
table(temp$M1, temp$state) #6 mismatches

temp <- markers[markers$M2Q > 1 & markers$state != "Unknown",]
table(temp$M2, temp$state) #10 mismatches

temp <- markers[markers$M3Q > 1 & markers$state != "Unknown",]
table(temp$M3, temp$state) #6 mismatches

#Changes over time
temp$period <- 0
temp$period <- ifelse(grepl("D",temp$ID),"Before", temp$period)
temp$period <- ifelse(grepl("E",temp$ID),"Before", temp$period)
temp$period <- ifelse(grepl("F",temp$ID),"After", temp$period)
temp$period <- ifelse(grepl("G",temp$ID),"After", temp$period)
temp$period <- ifelse(grepl("H",temp$ID),"After", temp$period)
temp$period <- ifelse(grepl("I",temp$ID),"After", temp$period)
temp$period <- ifelse(grepl("J",temp$ID),"After", temp$period)

table(temp$state, temp$period)


temp$event <- 0
temp$event <- ifelse(grepl("D",temp$ID),"D", temp$event)
temp$event <- ifelse(grepl("E",temp$ID),"E", temp$event)
temp$event <- ifelse(grepl("F",temp$ID),"F", temp$event)
temp$event <- ifelse(grepl("G",temp$ID),"G", temp$event)
temp$event <- ifelse(grepl("H",temp$ID),"H", temp$event)
temp$event <- ifelse(grepl("I",temp$ID),"I", temp$event)
temp$event <- ifelse(grepl("J",temp$ID),"J", temp$event)

table(temp$state, temp$period)
aggregate(temp$state, by=list(temp$period), FUN = length)
aggregate(temp$state, FUN = length)

table(temp$state)
?aggregate
#P20 RESISTANCE
#Load and format phenotype data
P20 <- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/time_series/IGV_files/IGVgroupsElocusAll.txt", header=TRUE)
P20 <- P20[,c(1,13)]
colnames(P20)[1] <- "ID"

markers <- merge(markers, P20, by = "ID")

temp <- markers[markers$phenotype == "RRS" | markers$phenotype == "RRR",]

table(temp$M1,temp$phenotype) # 6 mismatches
table(temp$M2,temp$phenotype) # 6 mismatches
table(temp$M3,temp$phenotype) # 5 mismatches
table(temp$state,temp$phenotype) #5 mismatches

#Manually import data for frequencies over time
freqs <- c(1.6,31.3,50.5,31.3,47.8,37.5)
haplo <- c("LSP 1","LSP 1","LSP 2","LSP 2","LSP 3","LSP 3")
time <- c(1,0,1,0,1,0)
df3 <- data.frame(haplo,time,freqs)

#Plot frequencies across time
ggplot(data=df3, aes(x=time, y=freqs,color = haplo))+
  geom_line()+
  geom_point(size=2)+
  theme_bw()+
  theme(panel.grid.major = element_blank() ,panel.grid.minor = element_blank())+
  ylab("Haplotype frequency")+ggtitle("Haplotype frequency changes")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_color_discrete(name = "Haplotype", labels = c("LSP1","LSP2","LSP3"))+
  xlim(-.10,1.10)+
  scale_x_continuous(breaks=c(0,1), labels = c("Pre-epidemic (n=56)", "Post-epidemic(n=186)"))

381x276
#Get number of samples in table
length(temp$phenotype)

#Pasteuria lineages
past <- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/time_series/IGV_files/PastDAPCgroups.txt")
markers <- merge(markers, past, by = "ID")

temp <- markers[markers$phenotype == "RRS" | markers$phenotype == "RRR",]
table(temp$group, temp$state)
table(temp$group, temp$M3)

table(temp$event)

#Put everything in a plot format
df2 <- temp
df2$X <- 0
df2$X <- ifelse(df2$state == "LSP1" | df2$state == "1/2" | df2$state == "1/3" ,"LSP1",df2$X)
df2$X <- ifelse(df2$state == "LSP2" | df2$state == "2/3" ,"LSP2",df2$X)
df2$X <- ifelse(df2$state == "LSP3" ,"LSP3",df2$X)

df2$Y <- 0
df2$Y <- ifelse(df2$state == "LSP1","LSP1",df2$Y)
df2$Y <- ifelse(df2$state == "LSP2" | df2$state == "1/2","LSP2",df2$Y)
df2$Y <- ifelse(df2$state == "LSP3" | df2$state == "1/3" | df2$state == "2/3","LSP3",df2$Y)

df2 <- df2[,c(1,13,15,16)]
#Import whole genomes from swisspond
genomes <- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/data/phenotypes/genome_assembly_phenotypes.txt", header=FALSE)
colnames(genomes) <- colnames(df2)
df2<- rbind(df2,genomes)
library(ggplot2)

p <- ggplot(data=df2, aes(x=X,y=Y))+
                         geom_jitter(height = 0.2,width = 0.2, aes(color=phenotype), size = 2, alpha = 0.7)+
                         theme_bw()+
                         theme(panel.grid.major = element_blank() ,panel.grid.minor = element_blank())+
                         xlab("LSP type")+ylab("LSP type")+ggtitle("P20 Susceptibility by LSP genotype")+
                         theme(plot.title = element_text(hjust = 0.5))+
                         scale_color_discrete(name = "Phenotype", labels = c("Resistant","Susceptible"))+
                         scale_x_discrete(labels = c("LSP1","LSP2","LSP3"))+
  theme(legend.position = c(0.8, 0.25),legend.background = element_rect(fill = "white", color = "black"))

ggsave("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/figures/geno_pheno_matrix.png",plot = p, width = 4, height = 3)

df$event <- 0

                       
                       