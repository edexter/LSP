### Minimal riparian plot

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