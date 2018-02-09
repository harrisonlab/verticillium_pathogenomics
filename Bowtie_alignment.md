##Alignment of raw reads vs the 12008 genome


### Vs the unmasked genome

```bash
Reference=$(ls repeat_masked/V.dahliae/12008/ncbi*/*_contigs_unmasked.fa)
for StrainPath in $(ls -d qc_dna/V.dahliae/paired/*/*)
do
echo $StrainPath
Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
echo "$Organism - $Strain"
F_Read=$(ls $StrainPath/F/*_trim.fq.gz | head -n1 | tail -n1)
R_Read=$(ls $StrainPath/R/*_trim.fq.gz | head -n1 | tail -n1)
echo $F_Read
echo $R_Read
OutDir=analysis/genome_alignment/bowtie/$Organism/$Strain/vs_12008unmasked
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment
qsub $ProgDir/bowtie/sub_bowtie.sh $Reference $F_Read $R_Read $OutDir $Strain
done
```

Identify read coverage over each bp

```bash
  for Bam in $(ls analysis/genome_alignment/bowtie/V.dahliae/*/vs_12008unmasked/*_aligned_sorted.bam); do
    Strain=$(echo $Bam | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Bam | rev | cut -f4 -d '/' | rev)
    echo "$Organism - $Strain"
    OutDir=$(dirname $Bam)
    # samtools depth -aa $Bam > $OutDir/${Organism}_${Strain}_vs_12008_depth.tsv
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment/coverage_analysis
    $ProgDir/cov_by_window.py --cov $OutDir/${Organism}_${Strain}_vs_12008_depth.tsv > $OutDir/${Organism}_${Strain}_vs_12008_depth_10kb.tsv
    sed -i "s/$/\t$Strain/g" $OutDir/${Organism}_${Strain}_vs_12008_depth_10kb.tsv
  done
  OutDir=analysis/genome_alignment/bowtie/grouped
  mkdir -p $OutDir
  cat analysis/genome_alignment/bowtie/*/*/*vs_12008unmasked/*_*_vs_12008_depth_10kb.tsv > analysis/genome_alignment/bowtie/grouped/vs_12008_grouped_depth.tsv
```



```R
library(readr)
appended_df <- read_delim("~/Downloads/Vd/vs_12008_grouped_depth.tsv",
     "\t", escape_double = FALSE, col_names = FALSE,
     trim_ws = TRUE)
colnames(appended_df) <- c("contig","position", "depth", "strain")
# appended_df$contig <- paste("Chr", appended_df$contig, sep = "")
appended_df$strain[appended_df$strain == "51"] <- "12251"
appended_df$strain[appended_df$strain == "53"] <- "12253"
appended_df$strain[appended_df$strain == "58"] <- "12158"
appended_df$strain[appended_df$strain == "61"] <- "12161"
appended_df$strain <- factor(appended_df$strain, levels = c("12008", "12251", "12253", "12158", "12161"))


# install.packages("ggplot2")
library(ggplot2)
require(scales)

for (i in 1:103){
  contig = paste("contig", i, sep = "_")
  p0 <- ggplot(data=appended_df[appended_df$contig == contig, ], aes(x=`position`, y=`depth`, group=1)) +
    geom_line() +
    labs(x = "Position (bp)", y = "Coverage") +
    scale_y_continuous(breaks=seq(0,150,25), limits=c(0,150)) +
    facet_wrap(~strain, nrow = 5, ncol = 1, strip.position = "left")
  outfile = paste("contig", i, "cov.tiff", sep = "_")
  ggsave(outfile , plot = p0, device = 'tiff', path = NULL,
    scale = 1, width = 250, height = 100, units = 'mm',
    dpi = 150, limitsize = TRUE)
  }
```




### Vs the hardmasked genome

```bash
  Reference=$(ls repeat_masked/V.dahliae/12008/ncbi_filtered_contigs_repmask/12008_contigs_hardmasked_repeatmasker_TPSI_appended.fa)
  for StrainPath in $(ls -d qc_dna/paired/V.dahliae/*)
  do
    echo $StrainPath
    Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
    Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
    echo "$Organism - $Strain"
    F_Read=$(ls $StrainPath/F/*_trim.fq.gz | head -n1 | tail -n1)
    R_Read=$(ls $StrainPath/R/*_trim.fq.gz | head -n1 | tail -n1)
    echo $F_Read
    echo $R_Read
    OutDir=analysis/genome_alignment/bwa/$Organism/$Strain/vs_12008_hardmasked
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment
    Prefix="${Strain}_vs_12008_hardmasked"
    # qsub $ProgDir/bowtie/sub_bowtie.sh $Reference $F_Read $R_Read $OutDir $Strain
    qsub $ProgDir/bwa/sub_bwa.sh $Prefix $Reference $F_Read $R_Read $OutDir
  done
```


Identify read coverage over each bp

```bash
  for Bam in $(ls analysis/genome_alignment/bowtie/V.dahliae/*/vs_12008_hardmasked/*_aligned_sorted.bam | grep -v '/12008/'); do
    Strain=$(echo $Bam | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Bam | rev | cut -f4 -d '/' | rev)
    echo "$Organism - $Strain"
    OutDir=$(dirname $Bam)
    # samtools depth -aa $Bam > $OutDir/${Organism}_${Strain}_vs_12008_depth.tsv
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment/coverage_analysis
    $ProgDir/cov_by_window.py --cov $OutDir/${Organism}_${Strain}_vs_12008_depth.tsv > $OutDir/${Organism}_${Strain}_vs_JR2_depth_10kb.tsv
    sed -i "s/$/\t$Strain/g" $OutDir/${Organism}_${Strain}_vs_12008_depth_10kb.tsv
  done
  OutDir=analysis/genome_alignment/bowtie/grouped
  mkdir -p $OutDir
  cat analysis/genome_alignment/bowtie/*/*/*vs_12008_hardmasked/*_*_vs_12008_depth_10kb.tsv > analysis/genome_alignment/bowtie/grouped/vs_12008_grouped_depth.tsv
  # paste analysis/genome_alignment/bowtie/V.dahliae/*/vs_12008unmasked/*_vs_12008_depth.tsv | cut -f 1,2,3,6,9,12,15 > $OutDir/vs_12008_grouped_depth.tsv
```


```bash
  for Bam in $(ls analysis/genome_alignment/bowtie/*/*/*vs_JR2/*_sorted.bam); do
    Strain=$(echo $Bam | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Bam | rev | cut -f4 -d '/' | rev)
    echo "$Organism - $Strain"
    OutDir=$(dirname $Bam)
    # samtools depth -aa $Bam > $OutDir/${Organism}_${Strain}_vs_JR2_depth.tsv
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment/coverage_analysis
    $ProgDir/cov_by_window.py --cov $OutDir/${Organism}_${Strain}_vs_JR2_depth.tsv > $OutDir/${Organism}_${Strain}_vs_JR2_depth_10kb.tsv
    sed -i "s/$/\t$Strain/g" $OutDir/${Organism}_${Strain}_vs_JR2_depth_10kb.tsv
  done
  cat analysis/genome_alignment/bowtie/*/*/*vs_JR2/*_*_vs_JR2_depth_10kb.tsv > analysis/genome_alignment/bowtie/grouped/vs_JR2_grouped_depth.tsv
```



```R
library(readr)


cov_depth <- read_delim("~/Downloads/Vd/tmp3.txt",
     "\t", escape_double = FALSE, col_names = FALSE,
     trim_ws = TRUE)
colnames(cov_depth) <- c("contig","position", "12008", "51", "53", "58", "61")


summary(cov_depth)

is.num <- sapply(cov_depth, is.numeric)
cov_depth[is.num] <- lapply(cov_depth[is.num], round, 0)

#install.packages("reshape")
library(reshape)
mdata <- reshape(cov_depth, direction='long', varying=c('12008', '51', '53', '58', '61'), v.names=c('strain', 'position'))

df_12008 <- cov_depth[,-7][, -6][, -5][, -4]
colnames(df_12008) <- c("contig","position", "depth")
df_12008$Strain <- "12008"

df_51 <- cov_depth[,-7][, -6][, -5][, -3]
colnames(df_51) <- c("contig","position", "depth")
df_51$Strain <- "51"

df_53 <- cov_depth[,-7][, -6][, -4][, -3]
colnames(df_53) <- c("contig","position", "depth")
df_53$Strain <- "53"

df_58 <- cov_depth[,-7][, -5][, -4][, -3]
colnames(df_58) <- c("contig","position", "depth")
df_58$Strain <- "58"

df_61 <- cov_depth[,-6][, -5][, -4][, -3]
colnames(df_61) <- c("contig","position", "depth")
df_61$Strain <- "61"

appended_df <- data.frame()
appended_df <- rbind(df_12008, df_51, df_53, df_58, df_61)

install.packages("ggplot2")
library(ggplot2)
require(scales)

p0 <- ggplot(data=appended_df[appended_df$contig == "contig_1", ], aes(x=`position`, y=`depth`, group=1)) +
  geom_line() +
  labs(x = "Position", y = "Coverage") +
  scale_y_continuous(breaks=seq(0,200,25), limits=c(0,200)) +
  facet_wrap(~Strain, nrow = 5, ncol = 1, strip.position = "left")

for (i in 1:103){
  contig = paste("contig_", i, sep = "")
  p0 <- ggplot(data=appended_df[appended_df$contig == contig, ], aes(x=`position`, y=`depth`, group=1)) +
    geom_line() +
    labs(x = "Position", y = "Coverage") +
    scale_y_continuous(breaks=seq(0,150,25), limits=c(0,150)) +
    facet_wrap(~Strain, nrow = 5, ncol = 1, strip.position = "left")
    outfile = paste(contig, "_cov.tiff", sep = "")
    ggsave(outfile , plot = p0, device = 'tiff', path = NULL,
      scale = 1, width = 100, height = 100, units = 'mm',
      dpi = 150, limitsize = TRUE)
  }

```



*************************************************************************************
Investigate the inversion event in 12008 genome. Firstly, check the assembly status of 12008 genome.

# Vs JR2

##Alignment of raw reads vs the JR2 genome

```bash
Reference=$(ls assembly/merged_canu_spades/V.dahliae/JR2/ensembl/Verticillium_dahliaejr2.GCA_000400815.2.dna.toplevel.fa)
for StrainPath in $(ls -d qc_dna/paired/V.dahliae/*/* | grep -e '12008' -e '51')
do
echo $StrainPath
Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
echo "$Organism - $Strain"
F_Read=$(ls $StrainPath/F/*_trim.fq.gz | head -n1 | tail -n1)
R_Read=$(ls $StrainPath/R/*_trim.fq.gz | head -n1 | tail -n1)
echo $F_Read
echo $R_Read
OutDir=analysis/genome_alignment/bowtie/$Organism/$Strain/vs_JR2
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment
qsub $ProgDir/bowtie/sub_bowtie.sh $Reference $F_Read $R_Read $OutDir $Strain
done
```

## Calculate coverage over the Ave1 region:

```bash
for Bam in $(ls analysis/genome_alignment/bowtie/*/*/*vs_JR2/*_sorted.bam); do
  Strain=$(echo $Bam | rev | cut -f3 -d '/' | rev)
  Organism=$(echo $Bam | rev | cut -f4 -d '/' | rev)
  echo "$Organism - $Strain"
  OutDir=$(dirname $Bam)
  Region="5:500000-1100000" # <chr:from-to>
  samtools depth -aa $Bam -r $Region > $OutDir/${Organism}_${Strain}_vs_JR2_depth.tsv
done
```
Download the files to my local machine and work in R studio:

```bash
scp cluster:/home/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics/analysis/genome_alignment/bowtie/V.dahliae/*/*vs_JR2/V.dahliae_*_vs_JR2_depth.tsv .
```

```R

library(readr)
V_dahliae_12008_vs_JR2_depth <- read_delim("~/Downloads/Vd/V.dahliae_12008_vs_JR2_depth.tsv",
     "\t", escape_double = FALSE, col_names = FALSE,
     trim_ws = TRUE)
colnames(V_dahliae_12008_vs_JR2_depth) <- c("contig","position", "depth")
V_dahliae_12008_vs_JR2_depth[,"Strain"]  <- "12008"
View(V_dahliae_12008_vs_JR2_depth)

V_dahliae_51_vs_JR2_depth <- read_delim("~/Downloads/Vd/V.dahliae_51_vs_JR2_depth.tsv",
     "\t", escape_double = FALSE, col_names = FALSE,
     trim_ws = TRUE)
colnames(V_dahliae_51_vs_JR2_depth) <- c("contig","position", "depth")
V_dahliae_51_vs_JR2_depth[,"Strain"]  <- "12251"

V_dahliae_53_vs_JR2_depth <- read_delim("~/Downloads/Vd/V.dahliae_53_vs_JR2_depth.tsv",
    "\t", escape_double = FALSE, col_names = FALSE,
    trim_ws = TRUE)
colnames(V_dahliae_53_vs_JR2_depth) <- c("contig","position", "depth")
V_dahliae_53_vs_JR2_depth[,"Strain"]  <- "12253"

V_dahliae_58_vs_JR2_depth <- read_delim("~/Downloads/Vd/V.dahliae_58_vs_JR2_depth.tsv",
   "\t", escape_double = FALSE, col_names = FALSE,
   trim_ws = TRUE)
colnames(V_dahliae_58_vs_JR2_depth) <- c("contig","position", "depth")
V_dahliae_58_vs_JR2_depth[,"Strain"]  <- "12158"

V_dahliae_61_vs_JR2_depth <- read_delim("~/Downloads/Vd/V.dahliae_61_vs_JR2_depth.tsv",
    "\t", escape_double = FALSE, col_names = FALSE,
    trim_ws = TRUE)
colnames(V_dahliae_61_vs_JR2_depth) <- c("contig","position", "depth")
V_dahliae_61_vs_JR2_depth[,"Strain"]  <- "12161"

appended = data.frame()
appended <- rbind(V_dahliae_12008_vs_JR2_depth, V_dahliae_51_vs_JR2_depth, V_dahliae_53_vs_JR2_depth, V_dahliae_58_vs_JR2_depth, V_dahliae_61_vs_JR2_depth)

appended<-appended[!(appended$position=="500000"),]
appended_mean = appended[seq(1, nrow(appended), 50), ]
appended_mean$depth <- colMeans(matrix(appended$depth, nrow=50))
appended_mean$Strain <- factor(appended_mean$Strain, levels = c("12008", "12251", "12253", "12158", "12161"))

require(scales)
library(ggplot2)

p0 <- ggplot(data=appended_mean, aes(x=`position`, y=`depth`, group=1)) +
  geom_line() +
  labs(x = "Position (bp)", y = "Coverage") +
  scale_x_continuous(breaks=seq(500000,1100000,100000), limits=c(500000,1100000)) +
  scale_y_continuous(breaks=seq(0,200,25), limits=c(0,200)) +
  facet_wrap(~Strain, nrow = 5, ncol = 1, strip.position = "left")
```


Identify read coverage over each bp

```bash
  for Bam in $(ls analysis/genome_alignment/bowtie/*/*/*vs_JR2/*_sorted.bam); do
    Strain=$(echo $Bam | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Bam | rev | cut -f4 -d '/' | rev)
    echo "$Organism - $Strain"
    OutDir=$(dirname $Bam)
    # samtools depth -aa $Bam > $OutDir/${Organism}_${Strain}_vs_JR2_depth.tsv
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment/coverage_analysis
    $ProgDir/cov_by_window.py --cov $OutDir/${Organism}_${Strain}_vs_JR2_depth.tsv > $OutDir/${Organism}_${Strain}_vs_JR2_depth_10kb.tsv
    sed -i "s/$/\t$Strain/g" $OutDir/${Organism}_${Strain}_vs_JR2_depth_10kb.tsv
  done
  cat analysis/genome_alignment/bowtie/*/*/*vs_JR2/*_*_vs_JR2_depth_10kb.tsv > analysis/genome_alignment/bowtie/grouped/vs_JR2_grouped_depth.tsv
```


```R
library(readr)
appended_df <- read_delim("~/Downloads/Vd/vs_JR2_grouped_depth.tsv",
     "\t", escape_double = FALSE, col_names = FALSE,
     trim_ws = TRUE)
colnames(appended_df) <- c("contig","position", "depth", "strain")
appended_df$contig <- paste("Chr", appended_df$contig, sep = "")
appended_df$strain[appended_df$strain == "51"] <- "12251"
appended_df$strain[appended_df$strain == "53"] <- "12253"
appended_df$strain[appended_df$strain == "58"] <- "12158"
appended_df$strain[appended_df$strain == "61"] <- "12161"
appended_df$strain <- factor(appended_df$strain, levels = c("12008", "12251", "12253", "12158", "12161"))


# install.packages("ggplot2")
library(ggplot2)
require(scales)

for (i in 1:8){
  contig = paste("Chr", i, sep = "")
  p0 <- ggplot(data=appended_df[appended_df$contig == contig, ], aes(x=`position`, y=`depth`, group=1)) +
    geom_line() +
    labs(x = "Position", y = "Coverage") +
    scale_y_continuous(breaks=seq(0,150,25), limits=c(0,150)) +
    facet_wrap(~strain, nrow = 5, ncol = 1, strip.position = "left")
  outfile = paste("chromosome", i, "cov.tiff", sep = "_")
  ggsave(outfile , plot = p0, device = 'tiff', path = NULL,
    scale = 1, width = 250, height = 100, units = 'mm',
    dpi = 150, limitsize = TRUE)
  }
```



# vs 12008 hardmasked genome
<!--
##Alignment of 12008 Miseq reads vs the 12008 genome to see the assembly status.

```bash
Reference=$(ls repeat_masked/V.dahliae/12008/ncbi*/*_contigs_unmasked.fa)
for StrainPath in $(ls -d qc_dna/paired/V.dahliae/12008)
do
echo $StrainPath
Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
echo "$Organism - $Strain"
F_Read=$(ls $StrainPath/F/*_trim.fq.gz | head -n1 | tail -n1)
R_Read=$(ls $StrainPath/R/*_trim.fq.gz | head -n1 | tail -n1)
echo $F_Read
echo $R_Read
OutDir=analysis/genome_alignment/bowtie/$Organism/$Strain/vs_12008unmasked
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment
qsub $ProgDir/bowtie/sub_bowtie.sh $Reference $F_Read $R_Read $OutDir $Strain
done
```

##Alignment of 12008 PacBio reads vs the 12008 genome to see the assembly status.
```bash
  Assembly=repeat_masked/V.dahliae/12008/ncbi*/*_contigs_unmasked.fa
  Reads=raw_dna/pacbio/V.dahliae/12008/extracted/concatenated_pacbio.fastq
  OutDir=analysis/genome_alignment/bwa/Verticillium/12008/PacReads_vs_12008merge
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment/bwa
  qsub $ProgDir/sub_bwa_pacbio.sh $Assembly $Reads $OutDir

```

##Alignment of 12008 PacBio reads vs the JR2 genome to see the assembly status.

```bash
Assembly=assembly/merged_canu_spades/V.dahliae/JR2/ensembl/Verticillium_dahliaejr2.GCA_000400815.2.dna.toplevel.fa
Reads=raw_dna/pacbio/V.dahliae/12008/extracted/concatenated_pacbio.fastq
OutDir=analysis/genome_alignment/bwa/Verticillium/12008/PacReads_vs_JR2
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment/bwa
qsub $ProgDir/sub_bwa_pacbio.sh $Assembly $Reads $OutDir
```

##Alignment of 12008 Miseq reads vs the JR2 genome to see the assembly status.

```bash
Reference=$(ls assembly/merged_canu_spades/V.dahliae/JR2/ensembl/Verticillium_dahliaejr2.GCA_000400815.2.dna.toplevel.fa)
for StrainPath in $(ls -d qc_dna/paired/V.dahliae/12008)
do
echo $StrainPath
Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
echo "$Organism - $Strain"
F_Read=$(ls $StrainPath/F/*_trim.fq.gz | head -n1 | tail -n1)
R_Read=$(ls $StrainPath/R/*_trim.fq.gz | head -n1 | tail -n1)
echo $F_Read
echo $R_Read
OutDir=analysis/genome_alignment/bowtie/$Organism/$Strain/Miseq_vs_JR2
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment
qsub $ProgDir/bowtie/sub_bowtie.sh $Reference $F_Read $R_Read $OutDir $Strain
done
```

##Alignment of 12008 Miseq reads vs the 12008 canu genome to revise the assembly status.

```bash
Reference=$(ls assembly/canu/V.dahliae/12008/polished/pilon.fasta)
for StrainPath in $(ls -d qc_dna/paired/V.dahliae/12008)
do
echo $StrainPath
Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
echo "$Organism - $Strain"
F_Read=$(ls $StrainPath/F/*_trim.fq.gz | head -n1 | tail -n1)
R_Read=$(ls $StrainPath/R/*_trim.fq.gz | head -n1 | tail -n1)
echo $F_Read
echo $R_Read
OutDir=analysis/genome_alignment/bowtie/$Organism/$Strain/Miseq_vs_12008canu
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment
qsub $ProgDir/bowtie/sub_bowtie.sh $Reference $F_Read $R_Read $OutDir $Strain
done
```

##Alignment of 12008 Miseq reads vs the 12008 Spades genome to revise the assembly status.

```bash
Reference=$(ls assembly/spades_pacbio/V.dahliae/12008/filtered_contigs/contigs_min_500bp.fasta)
for StrainPath in $(ls -d qc_dna/paired/V.dahliae/12008)
do
echo $StrainPath
Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
echo "$Organism - $Strain"
F_Read=$(ls $StrainPath/F/*_trim.fq.gz | head -n1 | tail -n1)
R_Read=$(ls $StrainPath/R/*_trim.fq.gz | head -n1 | tail -n1)
echo $F_Read
echo $R_Read
OutDir=analysis/genome_alignment/bowtie/$Organism/$Strain/Miseq_vs_12008spades
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment
qsub $ProgDir/bowtie/sub_bowtie.sh $Reference $F_Read $R_Read $OutDir $Strain
done
``` -->

```bash


```
