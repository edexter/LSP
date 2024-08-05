# Genespace riparian plot chromosome 6

````
R

library(GENESPACE)

load('/scicore/home/ebertd/dexter0000/LSP/genespace/working_dir/results/gsParams.rda',verbose = TRUE)

CH_434_inb3_a_1<-data.frame(
  genome = c("CH_434_inb3_a_1"), 
  chr = c(12,37,11,17))

NO_V_7<-data.frame(
  genome = c("NO_V_7"), 
  chr = c(71,2,77,69,11))
  
FI_SK_58_2_18_4<-data.frame(
  genome = c("FI_SK_58_2_18_4"), 
  chr = c(4,8,25,42,13))
  
Xinb3<-data.frame(
  genome = c("Xinb3"), 
  chr = c(1,108,49,80,61))

RU_RM1_2<-data.frame(
  genome = c("RU_RM1_2"), 
  chr = c(1,70,73))
  
CN_W1_1<-data.frame(
  genome = c("CN_W1_1"), 
  chr = c(35,21,49,16,69,0,71,13,9))
 
t3_12_3_1i_12<-data.frame(
  genome = c("t3_12_3_1i_12"), 
  chr = c(33,20,17,75,157,10,137,76))

CH_t4_12_3_3<-data.frame(
  genome = c("CH_t4_12_3_3"), 
  chr = c(64,8,50,3,129))
  
CH_H_2299<-data.frame(
genome = c("CH_H_2299"), 
chr = c(85,97,29,24,38,25))

CH_H_2015_59<-data.frame(
  genome = c("CH_H_2015_59"), 
  chr = c(41,86,147,92,14,25,31))

US_SP_221_1<-data.frame(
  genome = c("US_SP_221_1"), 
  chr = c(37))

D_similis_IL_SIM_A20<-data.frame(
  genome = c("D_similis_IL_SIM_A20"), 
  chr = c(134,256,338,104,111,81,43,64,93,69,17))

CH_H_2015_49<-data.frame(
  genome = c("CH_H_2015_49"), 
  chr = c(76,35))

t2_17_3_4i_13<-data.frame(
  genome = c("t2_17_3_4i_13"),
    chr = c(191,46,161))

t1_10_3_2<-data.frame(
  genome = c("t1_10_3_2"), 
  chr = c(161,23,35,261,30,248,47,56,107,180,482,5,218,64,159,156,155))

NCBI_scaffold<-data.frame(
  genome = c("NCBI_scaffold"), 
  chr = c(361461,100000111))

CH14_scaffold<-data.frame(
  genome = c("CH14_scaffold"), 
  chr = c(9))

invchr <- rbind(CH_434_inb3_a_1, NO_V_7, FI_SK_58_2_18_4, Xinb3, RU_RM1_2, CN_W1_1, t3_12_3_1i_12, CH_H_2299, CH_H_2015_59, US_SP_221_1, D_similis_IL_SIM_A20, t2_17_3_4i_13, t1_10_3_2, CH_H_2015_49,NCBI_scaffold,CH14_scaffold)

genomeIds <- c("D_similis_IL_SIM_A20","US_SP_221_1","RU_RM1_2","CN_W1_1","CH_H_2299","CH_434_inb3_a_1","FI_SK_58_2_18_4","CH_H_2015_49","t3_12_3_1i_12","CH_H_2015_59","t2_17_3_4i_13")

roi <- data.frame(
  genome = c("t2_17_3_4i_13"), 
  chr = c("151","131"))

pdf(file = "/scicore/home/ebertd/dexter0000/LSP/genespace/working_dir/CHR6b.pdf",
    width = 8,
    height = 9)
    
ripDat <- plot_riparian(
  gsParam, 
  refGenome = "t2_17_3_4i_13",
  useOrder = FALSE, 
  useRegions = FALSE,
  reorderBySynteny = TRUE,
  #addThemes = ggthemes,
  syntenyWeight = 1,
  #customRefChrOrder = c("15_1","13_1"),
  #backgroundColor = NULL,
  xlabel = sprintf("D. magna chromosome 6"),
  #inversionColor = "white",
  genomeIDs = genomeIds,
  highlightBed = roi, 
  backgroundColor = NULL,
  invertTheseChrs = invchr,
  minChrLen2plot = 50000
)
  
  dev.off()
````

