# Create initial riparian plot with all genomes

````R
srun --nodes=1 --cpus-per-task=2 --mem=16G --pty bash

cd /scicore/home/ebertd/dexter0000/LSP

conda activate orthofinder

module load R

R

library(GENESPACE)

#Load genespace results
load('/scicore/home/ebertd/dexter0000/LSP/genespace/working_dir/results/gsParams.rda',verbose = TRUE)
        
################################################################
#Manual curation of genomes
################################################################
CH_434_inb3_a_1<-data.frame(
  genome = c("CH_434_inb3_a_1"), 
  chr = c(12,37,5))

NO_V_7<-data.frame(
  genome = c("NO_V_7"), 
  chr = c(77,69,11))
  
FI_SK_58_2_18_4<-data.frame(
  genome = c("FI_SK_58_2_18_4"), 
  chr = c(8,39,31,4,42,25))
  
Xinb3<-data.frame(
  genome = c("Xinb3"), 
  chr = c(108,49,80,61))

RU_RM1_2<-data.frame(
  genome = c("RU_RM1_2"), 
  chr = c(1,70,6,114,73))
  
CN_W1_1<-data.frame(
  genome = c("CN_W1_1"), 
  chr = c(35,21,49,16,69,0))
 
t3_12_3_1i_12<-data.frame(
  genome = c("t3_12_3_1i_12"), 
  chr = c(75,33,157,20,17))

CH_t4_12_3_3<-data.frame(
  genome = c("CH_t4_12_3_3"), 
  chr = c(26,76,29,3,43))
  
CH_H_2299<-data.frame(
genome = c("CH_H_2299"), 
chr = c(77,85,97))

CH_H_2015_59<-data.frame(
  genome = c("CH_H_2015_59"), 
  chr = c(41,86,147,92,14,25,108,31))

US_SP_221_1<-data.frame(
  genome = c("US_SP_221_1"), 
  chr = c(37))

D_similis_IL_SIM_A20<-data.frame(
  genome = c("D_similis_IL_SIM_A20"), 
  chr = c(134,256,448,64,81,43,31,338,104,111))

CH_H_2015_49<-data.frame(
  genome = c("CH_H_2015_49"), 
  chr = c(68,76,35,151))

t2_17_3_4i_13<-data.frame(
  genome = c("t2_17_3_4i_13"),
    chr = c(46))

t1_10_3_2<-data.frame(
  genome = c("t1_10_3_2"), 
  chr = c(161,23,208,82,662,35,56,98))

CH_2016_H_34<-data.frame(
  genome = c("CH_2016_H_34"), 
  chr = c(55,164,21,74,105))

US_D_3<-data.frame(
  genome = c("US_D_3"), 
  chr = c(126,7,98,27,38,106,77,9))

NCBI_scaffold<-data.frame(
  genome = c("NCBI_scaffold"), 
  chr = c(361461))

ET_C_1<-data.frame(
  genome = c("ET_C_1"), 
  chr = c(28,64,2))

DZ_JV_2<-data.frame(
  genome = c("DZ_JV_2"), 
  chr = c(26,189,73,116,75,103,19))

IL_TY_10<-data.frame(
  genome = c("IL_TY_10"), 
  chr = c(25,56,15))

invchr <- rbind(CH_434_inb3_a_1, NO_V_7, FI_SK_58_2_18_4, Xinb3, RU_RM1_2, CN_W1_1, t3_12_3_1i_12, t3_12_3_1i_12, CH_t4_12_3_3, CH_H_2299, CH_H_2015_59, US_SP_221_1, D_similis_IL_SIM_A20, t2_17_3_4i_13, t1_10_3_2, CH_H_2015_49,NCBI_scaffold,CH_2016_H_34,US_D_3,ET_C_1,DZ_JV_2,IL_TY_10)

genomeIds <- c("D_similis_IL_SIM_A20",             
               "US_SP_221_1",
               "RU_RM1_2",
               "CN_W1_1",               
               "CH_H_2299",               
               "CH_434_inb3_a_1",
               "FI_SK_58_2_18_4",               
               "CH_H_2015_59",              
               "t3_12_3_1i_12",
               "CH_H_2015_49",
               "t2_17_3_4i_13"
              )

##############################################################################
#plots
##############################################################################
roi <- data.frame(
  genome = c("t2_17_3_4i_13"), 
  chr = c("191","161","31","351"))

# export plot
png(file = "/scicore/home/ebertd/dexter0000/LSP/genespace/working_dir/riparian.png",
    width = 8,
    height = 9, units="in", res = 900)
    
ripDat <- plot_riparian(
  gsParam, 
  refGenome = "t2_17_3_4i_13",
  useOrder = FALSE, 
  useRegions = FALSE,
  reorderBySynteny = TRUE,
  syntenyWeight = 1,
  xlabel = sprintf("D. magna chromosome 5"),
  genomeIDs = genomeIds,
  highlightBed = roi, 
  backgroundColor = NULL,
  invertTheseChrs = invchr,
  minChrLen2plot = 50000,
  forceRecalcBlocks = FALSE,
  scalePlotWidth = 1,
  gapProp = 0.002,
  braidAlpha = 0.8,
  customRefChrOrder=c("191","161","31","351")
)
  
  dev.off()

````



# Pairwise just three

````
genomeIds <- c("t1_10_3_2",             
               "CH_434_inb3_a_1",
               "t2_17_3_4i_13"
              )

##############################################################################
#plots
##############################################################################
roi <- data.frame(
  genome = c("t2_17_3_4i_13"), 
  chr = c("191","161","31","351"))

# export plot
pdf(file = "/scicore/home/ebertd/dexter0000/LSP/genespace/working_dir/riparian_swisspond.pdf",
    width = 7,
    height = 9)
    
ripDat <- plot_riparian(
  gsParam, 
  refGenome = "t2_17_3_4i_13",
  useOrder = FALSE, 
  useRegions = FALSE,
  reorderBySynteny = TRUE,
  syntenyWeight = 1,
  xlabel = sprintf("D. magna chromosome 5"),
  genomeIDs = genomeIds,
  highlightBed = roi, 
  backgroundColor = NULL,
  invertTheseChrs = invchr,
  minChrLen2plot = 50000,
  forceRecalcBlocks = FALSE,
  scalePlotWidth = 1,
  gapProp = 0.002,
  braidAlpha = 0.8,
  customRefChrOrder=c("191","161","31","351")
)
  
  dev.off()
````









