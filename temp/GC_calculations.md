# Calculate GC % across the genomes and contigs of interest

This is a quick calculation and can be performed interactively as needed using the script "nucleotide_biases.pl", which is given below.

````bash
srun --nodes=1 --cpus-per-task=32 --mem=4G --pty bash

cd /scicore/home/ebertd/dexter0000/LSP

# Haplotype 1 - chromsome 5 right arm contigs
mkdir GC/CH14_chr5.ptg000016l_1 GC/CH14_chr5.ptg000019l_1
perl scripts/nucleotide_biases.pl \
	-f haplotypes/CH14_chr5.fa \
	-w 50000 \
	-o GC \
	-s 1000

# Haplotype 2 - chromsome 5 right arm contigs
mkdir GC/CH434_chr5.000012F GC/CH434_chr5.000037F
perl scripts/nucleotide_biases.pl \
	-f haplotypes/CH434_chr5.fa \
	-w 50000 \
	-o GC \
	-s 1000

# Haplotype 3 - chromsome 5 right arm contigs
mkdir GC/T1_chr5.ptg000005l GC/T1_chr5.ptg000023l GC/T1_chr5.ptg000035l GC/T1_chr5.ptg000082l GC/T1_chr5.ptg000107l GC/T1_chr5.ptg000180l
perl scripts/nucleotide_biases.pl \
	-f haplotypes/T1_chr5.fa \
	-w 50000 \
	-o GC \
	-s 1000
    
````



### nucleotide_biases.pl

Here is the script (from the Pombert lab) that is used to calculate GC biases.

````perl
#!/usr/bin/perl
## Pombert Lab, 2022
my $name = 'nucleotide_biases.pl';
my $version = '0.2';
my $updated = '2022-06-22';

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Basename;
use File::Path qw(make_path);

my $usage = <<"OPTIONS";
NAME		$name
VERSION		$version
UPDATED		$updated
SYNOPSIS	Generates tab-delimited sliding windows of GC, AT, purine and pyrimidine
		 distributions for easy plotting with MS Excel or other tool.

COMMAND		$name \\
		  -f *.fasta \\
		  -o output_directory \\
		  -w 1000 \\
		  -s 500

-f (--fasta)	Fasta file(s) to process
-o (--outdir)	Output directory [Default: ntBiases]
-w (--winsize)	Sliding window size [Default: 1000]
-s (--step)		Sliding window step [Default: 500]
OPTIONS
die "\n$usage\n" unless @ARGV;

my @fasta;
my $outdir = 'ntBiases';
my $winsize = 1000;
my $step = 500;
GetOptions(
	'f|fasta=s@{1,}' => \@fasta,
	'o|outdir=s' => \$outdir,
	'w|winsize=i' => \$winsize,
	's|step=i' => \$step
);

### Check if output directory / subdirs can be created
$outdir =~ s/\/$//;
unless (-d $outdir) {
	make_path( $outdir, { mode => 0755 } )  or die "Can't create $outdir: $!\n";
}

### Iterating through FASTA file(s)
while (my $fasta = shift@fasta){

	my ($basename) = fileparse($fasta);
	my ($fileprefix) = $basename =~ /(\S+)\.\w+$/;

	open FASTA, "<", $fasta or die "Can't open $fasta: $!\n";

	### Creating database of sequences (could be multifasta)
	my %sequences;
	my $seqname;

	while (my $line = <FASTA>){
		chomp $line;
		if ($line =~ /^>(\S+)/){
			$seqname = $1;
		}
		else {
			$sequences{$seqname} .= $line;
		}
	}

	### Iterating through each sequence in the FASTA file
	foreach my $sequence (sort (keys %sequences)){

		my $outfile = $outdir.'/'.$fileprefix.'.'.$sequence.'.tsv';
		open GC, ">", $outfile or die "Can't create $outfile: $!\n";

		print GC "# Location\t% GC\t% AT\t% Purines\t% Pyrimidines\t% GT\t% AC\n";

		### Sliding windows
		my $seq = $sequences{$sequence};
		my $csize = length $seq;
		my $x;
		for ($x = 0; $x <= ($csize - $step); $x += $step){
			my $subseq = substr($seq, $x, $winsize);
			my $gc = $subseq =~ tr/GgCc//;
			my $at = $subseq =~ tr/AaTt//;
			my $pur = $subseq =~ tr/GgAa//;
			my $pyr = $subseq =~ tr/CcTt//;
			my $gt = $subseq =~ tr/GgTt//;
			my $ac = $subseq =~ tr/AaCc//;
			$gc = ($gc)/$winsize * 100;
			$at = ($at)/$winsize * 100;
			$pur = ($pur)/$winsize * 100;
			$pyr = ($pyr)/$winsize * 100;
			$gt = ($gt)/$winsize * 100;
			$ac = ($ac)/$winsize * 100;
			print GC "$x\t$gc\t$at\t$pur\t$pyr\t$gt\t$ac\n";
		}

		### Working on leftover string < $winsize
		my $modulo = $csize % $winsize;
		my $subseqleft = substr($seq, -$modulo, $modulo);
		my $leftover_size = length $subseqleft;
		my $gc = $subseqleft =~ tr/GgCc//;
		my $at = $subseqleft =~ tr/AaTt//;
		my $pur = $subseqleft =~ tr/GgAa//;
		my $pyr = $subseqleft =~ tr/CcTt//;
		my $gt = $subseqleft =~ tr/GgTt//;
		my $ac = $subseqleft =~ tr/AaCc//;
		$gc = ($gc)/$leftover_size * 100;
		$at = ($at)/$winsize * 100;
		$pur = ($pur)/$winsize * 100;
		$pyr = ($pyr)/$winsize * 100;
		$gt = ($gt)/$winsize * 100;
		$ac = ($ac)/$winsize * 100;
		print GC "$x\t$gc\t$at\t$pur\t$pyr\t$gt\t$ac\n";

		close GC;

	}

	close FASTA;

}
````



# R code for plotting

Here is the R code that can be run interactively to make the GC plots/

````R
################################################################################
# Load required packages
library(ggplot2)
library(gridExtra)

################################################################################
# Define global variables
scaleFactor <- 10^6
windowSize <- 50000

################################################################################
# GC plot for entire chromosome CH14

# Convert position in contig to position in chromosome
c19 <- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/data/GC/CH14/CH14_chr5.ptg000019l_1/rc.tsv")
cutoff <- max(c19$X..Location) - windowSize
c19 <- c19[c19$X..Location < cutoff,]
maxPos <- max(c19$X..Location)

# Convert position in contig to position in chromosome
c16 <- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/data/GC/CH14/CH14_chr5.ptg000016l_1/rc.tsv")
cutoff <- max(c16$X..Location) - windowSize
c16 <- c16[c16$X..Location < cutoff,]
c16$X..Location <- c16$X..Location + maxPos
maxPos <- max(c16$X..Location)

# Convert position in contig to position in chromosome
c3 <- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/data/GC/CH14/CH14_chr5.ptg000003l_1.tsv")
cutoff <- max(c3$X..Location) - windowSize
c3 <- c3[c3$X..Location < cutoff,]
c3$X..Location <- c3$X..Location + maxPos
maxPos <- max(c3$X..Location)

# Convert position in contig to position in chromosome
c35 <- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/data/GC/CH14/CH14_chr5.ptg000035l_1.tsv")
cutoff <- max(c35$X..Location) - windowSize
c35 <- c35[c35$X..Location < cutoff,]
c35$X..Location <- c35$X..Location + maxPos
maxPos <- max(c35$X..Location)

CH14 <- rbind(c19,c16,c3, c35)

# Generate plot
p <- ggplot(CH14,aes(x=X..Location / scaleFactor, y=X..GC))+
  geom_line(size=0.2)+
  geom_hline(yintercept = 40, col = "red", lty =2)+
  labs(x="Position in chromosome (Mb)", y="GC %")+
  theme_bw()+
  scale_x_continuous(limits = c(0, maxPos/scaleFactor), expand = c(0, 0))

p

#Save plot
ggsave(plot = p, "C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/figures/GC_CH14.png", width = 12, height = 6, units = "cm", bg = "white")

################################################################################
#Plot GC for just contig containing the LSP across multiple genomes

# Load haplotype 1 GC values
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

# Load haplotype 2 GC values
df <- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/data/GC/CH14/CH434_chr5.000012F/rc.tsv")
cutoff <- max(df$X..Location) - windowSize
df <- df[df$X..Location < cutoff,]
pB <- ggplot(df,aes(x=X..Location/scaleFactor, y=X..GC))+
  geom_line()+
  geom_hline(yintercept = 40, col = "blue", lty =2)+
  labs(x="Position in contig (Mb)", y="GC %")+
  theme_bw()+
  xlim(0,9)+ylim(25,50)+
  annotate("rect", xmin=537051/scaleFactor, xmax=3078906/scaleFactor, ymin=-Inf, ymax=Inf, alpha=0.25, fill="gray")+
  labs(title = "Haplotype 2 (2.5 Mb)")

# Load haplotype 3 GC values
df <- read.delim("C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/data/GC/CH14/T1_chr5.ptg000023l/rc.tsv")
cutoff <- max(df$X..Location) - windowSize
df <- df[df$X..Location < cutoff,]
pC <- ggplot(df,aes(x=X..Location/scaleFactor, y=X..GC))+
  geom_line()+
  geom_hline(yintercept = 40, col = "blue", lty =2)+
  labs(x="Position in contig (Mb)", y="GC %")+
  theme_bw()+
  xlim(0,9)+ ylim(25,50)+
  annotate("rect", xmin=627218/scaleFactor, xmax=3693875/scaleFactor, ymin=-Inf, ymax=Inf, alpha=0.25, fill="gray")+
  labs(title = "Haplotype 3 (3.1 Mb)",x="")

# Plot all 3 haplotypes
pAll <- grid.arrange(pA, pB, pC, nrow = 3)
pAll

# Save plot
ggsave(plot = pAll, "C:/Users/ericd/Dropbox/Eric Work/Ebert lab/LSP/figures/GCplots.png", width = 16, height = 12, units = "cm", bg = "white")
````

