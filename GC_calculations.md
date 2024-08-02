# Calcuate GC % across the genomes and contigs of interest

````bash
srun --nodes=1 --cpus-per-task=32 --mem=4G --pty bash

cd /scicore/home/ebertd/dexter0000/LSP

mkdir GC/CH14_chr5.ptg000016l_1 GC/CH14_chr5.ptg000019l_1
perl scripts/nucleotide_biases.pl \
	-f haplotypes/CH14_chr5.fa \
	-w 50000 \
	-o GC \
	-s 1000

mkdir GC/CH434_chr5.000012F GC/CH434_chr5.000037F
perl scripts/nucleotide_biases.pl \
	-f haplotypes/CH434_chr5.fa \
	-w 50000 \
	-o GC \
	-s 1000

mkdir GC/T1_chr5.ptg000005l GC/T1_chr5.ptg000023l GC/T1_chr5.ptg000035l GC/T1_chr5.ptg000082l GC/T1_chr5.ptg000107l GC/T1_chr5.ptg000180l
perl scripts/nucleotide_biases.pl \
	-f haplotypes/T1_chr5.fa \
	-w 50000 \
	-o GC \
	-s 1000
    
````

# R code for plotting

````R
library(ggplot2)
library(gridExtra)

scaleFactor <- 10^6
windowSize <- 50000

# Entire chromosome CH14
c19 <- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/data/GC/CH14/CH14_chr5.ptg000019l_1/rc.tsv")
cutoff <- max(c19$X..Location) - windowSize
c19 <- c19[c19$X..Location < cutoff,]
maxPos <- max(c19$X..Location)

c16 <- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/data/GC/CH14/CH14_chr5.ptg000016l_1/rc.tsv")
cutoff <- max(c16$X..Location) - windowSize
c16 <- c16[c16$X..Location < cutoff,]
c16$X..Location <- c16$X..Location + maxPos
maxPos <- max(c16$X..Location)

c3 <- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/data/GC/CH14/CH14_chr5.ptg000003l_1.tsv")
cutoff <- max(c3$X..Location) - windowSize
c3 <- c3[c3$X..Location < cutoff,]
c3$X..Location <- c3$X..Location + maxPos
maxPos <- max(c3$X..Location)

c35 <- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/data/GC/CH14/CH14_chr5.ptg000035l_1.tsv")
cutoff <- max(c35$X..Location) - windowSize
c35 <- c35[c35$X..Location < cutoff,]
c35$X..Location <- c35$X..Location + maxPos
maxPos <- max(c35$X..Location)

CH14 <- rbind(c19,c16,c3, c35)

p <- ggplot(CH14,aes(x=X..Location / scaleFactor, y=X..GC))+
  geom_line(size=0.2)+
  geom_hline(yintercept = 40, col = "red", lty =2)+
  labs(x="Position in chromosome (Mb)", y="GC %")+
  theme_bw()+
  scale_x_continuous(limits = c(0, maxPos/scaleFactor), expand = c(0, 0))

p

ggsave(plot = p, "C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/figures/GC_CH14.png", width = 12, height = 6, units = "cm", bg = "white")

################################################################################
#Plot GC bias for just contig containing the LSP
################################################################################

df <- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/data/GC/CH14/CH14_chr5.ptg000003l_1.tsv") 
cutoff <- max(df$X..Location) - windowSize
df <- df[df$X..Location < cutoff,]
pA <- ggplot(df,aes(x=X..Location/scaleFactor, y=X..GC))+
  geom_line()+
  geom_hline(yintercept = 40, col = "blue", lty =2)+
  labs(x="Position in contig (Mb)", y="GC %")+
  theme_bw()+
  xlim(0,9)+ ylim(25,50)+
  annotate("rect", xmin=1534355/scaleFactor, xmax=7111581/scaleFactor, ymin=-Inf, ymax=Inf, alpha=0.25, fill="gray")+
  labs(title = "Haplotype 1 (5.6 Mb)", x="")

# T1
df <- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/data/GC/CH14/T1_chr5.ptg000023l/rc.tsv")
cutoff <- max(df$X..Location) - windowSize
df <- df[df$X..Location < cutoff,]
pB <- ggplot(df,aes(x=X..Location/scaleFactor, y=X..GC))+
  geom_line()+
  geom_hline(yintercept = 40, col = "blue", lty =2)+
  labs(x="Position in contig (Mb)", y="GC %")+
  theme_bw()+
  xlim(0,9)+ ylim(25,50)+
  annotate("rect", xmin=627218/scaleFactor, xmax=3693875/scaleFactor, ymin=-Inf, ymax=Inf, alpha=0.25, fill="gray")+
  labs(title = "Haplotype 3 (3.1 Mb)",x="")

# CH434
df <- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/data/GC/CH14/CH434_chr5.000012F/rc.tsv")
cutoff <- max(df$X..Location) - windowSize
df <- df[df$X..Location < cutoff,]
pC <- ggplot(df,aes(x=X..Location/scaleFactor, y=X..GC))+
  geom_line()+
  geom_hline(yintercept = 40, col = "blue", lty =2)+
  labs(x="Position in contig (Mb)", y="GC %")+
  theme_bw()+
  xlim(0,9)+ylim(25,50)+
  annotate("rect", xmin=537051/scaleFactor, xmax=3078906/scaleFactor, ymin=-Inf, ymax=Inf, alpha=0.25, fill="gray")+
  labs(title = "Haplotype 2 (2.5 Mb)")

pAll <- grid.arrange(pA, pC, pB, nrow = 3)
pAll
ggsave(plot = pAll, "C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/figures/GCplots.png", width = 16, height = 12, units = "cm", bg = "white")
````

