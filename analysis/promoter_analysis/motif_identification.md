# Motif identification Vd

Commands used to identify motifs in the promoters of key groups of Vd DEGs

Firstly, promoter sequences were extracted from the genome using
published gene models

```bash
# Assembly=$(ls ../clocks/public_genomes/JR2/Verticillium_dahliaejr2.GCA_000400815.2.dna.toplevel_parsed.fa)
Assembly=$(ls ../clocks/public_genomes/JR2/Verticillium_dahliaejr2.GCA_000400815.2.dna.toplevel.fa)
Gff=$(ls ../clocks/public_genomes/JR2/Verticillium_dahliaejr2.GCA_000400815.2.33.gff3)

OutDir=analysis/promoters/V.dahliae/JR2
mkdir -p $OutDir
ProgDir=/home/armita/git_repos/emr_repos/scripts/verticillium_pathogenomics/analysis/promoter_analysis
$ProgDir/extract_promoters_JR2.py --gff $Gff --fasta $Assembly --prefix $OutDir/JR2_promoters

ls $OutDir
```


## WC-1 characterised genes

Motif were identified in the promoters of key groups of DEGs.

DEGs were first selected by:
* eye from heatmaps of DEGs under different conditions
* based on cluster analysis.


### Selected genes:

Attempt to identify common elements in manually selected groups of similarly regulated transcription factors.

```bash
screen -a
qlogin


ProjDir=/home/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics
cd $ProjDir
ProjDir=$(ls -d /home/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics)
Promoters=$(ls $ProjDir/analysis/promoters/V.dahliae/JR2/JR2_promoters.fa)
OutDir=analysis/promoters/motifs/V.dahliae/JR2/WC-1/selected_genes/selected
mkdir -p $OutDir/meme

# Transcription factors each showing decresed expression in frq knock outs

GeneList="VDAG_JR2_Chr2g00510 VDAG_JR2_Chr2g08700 VDAG_JR2_Chr2g10990 VDAG_JR2_Chr5g01930 VDAG_JR2_Chr3g13030 VDAG_JR2_Chr2g10850 VDAG_JR2_Chr1g24030 VDAG_JR2_Chr3g06490 VDAG_JR2_Chr3g06100 VDAG_JR2_Chr2g03270 VDAG_JR2_Chr5g01890 VDAG_JR2_Chr8g07070 VDAG_JR2_Chr4g04620 VDAG_JR2_Chr8g07580 VDAG_JR2_Chr1g27370 VDAG_JR2_Chr1g03770 VDAG_JR2_Chr3g12060 VDAG_JR2_Chr2g00520 VDAG_JR2_Chr8g00200 VDAG_JR2_Chr1g26850 VDAG_JR2_Chr1g12350 VDAG_JR2_Chr2g10410 VDAG_JR2_Chr6g02190 VDAG_JR2_Chr7g07100 VDAG_JR2_Chr5g05720 VDAG_JR2_Chr2g01990 VDAG_JR2_Chr8g01630 VDAG_JR2_Chr5g00990 VDAG_JR2_Chr1g23950 VDAG_JR2_Chr2g07860 VDAG_JR2_Chr5g11380 VDAG_JR2_Chr2g07010 VDAG_JR2_Chr1g27300 VDAG_JR2_Chr5g05660 VDAG_JR2_Chr3g09510 VDAG_JR2_Chr2g00700 VDAG_JR2_Chr1g07150 VDAG_JR2_Chr3g08070 VDAG_JR2_Chr3g08730 VDAG_JR2_Chr2g04990
VDAG_JR2_Chr1g08930 VDAG_JR2_Chr1g08230 VDAG_JR2_Chr2g05130 VDAG_JR2_Chr3g04790 VDAG_JR2_Chr1g05210 VDAG_JR2_Chr7g10270 VDAG_JR2_Chr3g03860 VDAG_JR2_Chr8g11340 VDAG_JR2_Chr6g09990 VDAG_JR2_Chr1g09490 VDAG_JR2_Chr1g18820 VDAG_JR2_Chr2g11140 VDAG_JR2_Chr5g04520 VDAG_JR2_Chr6g07650 VDAG_JR2_Chr1g17310"
ListLen=$(echo $GeneList | grep -o 'g' | wc -l)
SubPromoters=$OutDir/promoters_${ListLen}.fa

for GeneID in $GeneList; do
cat $Promoters | sed -n "/^>${GeneID}_/,/^>/p" | grep -v "^$" | head -n -1
done > $SubPromoters

GeneList="VDAG_JR2_Chr3g10380 VDAG_JR2_Chr1g12450 VDAG_JR2_Chr1g15260 VDAG_JR2_Chr3g07640 VDAG_JR2_Chr2g01990 VDAG_JR2_Chr1g17210"


ListLen=$(echo $GeneList | grep -o 'g' | wc -l)
SubPromoters=$OutDir/promoters_${ListLen}.fa

for GeneID in $GeneList; do
cat $Promoters | sed -n "/^>${GeneID}_/,/^>/p" | grep -v "^$" | head -n -1
done > $SubPromoters


meme $SubPromoters -dna -mod anr -nmotifs 5 -minw 6 -maxw 12 -revcomp -evt 1.0e+005 -oc $OutDir/meme_${ListLen}

cat $OutDir/meme_${ListLen}/meme.txt | grep -C2 'regular expression'

#---
# Running MAST
#---

mast $OutDir/meme_${ListLen}/meme.xml $SubPromoters -oc $OutDir/mast_${ListLen}
ls $OutDir/mast_${ListLen}/mast.txt

#---
# Running FIMO
#---
Thresh=0.00006 #CASSIS (trained for secmet clusters)
# Thresh=1e-4 #Default
fimo -thresh $Thresh -oc $OutDir/fimo_${ListLen} $OutDir/meme_${ListLen}/meme.html $OutDir/promoters_50.fa
for Motif in $(cat $OutDir/fimo_${ListLen}/fimo.tsv | grep 'TSS' | cut -f1 | sort | uniq); do
  echo $Motif
  echo "Promoters containing motif:"
  cat $OutDir/fimo_${ListLen}/fimo.tsv | grep "^$Motif" | cut -f3 | sort | uniq | wc -l
  echo "Total number found in these promoters:"
  cat $OutDir/fimo_${ListLen}/fimo.tsv | grep "^$Motif" | cut -f3 | sort | wc -l
  echo "this covers the following genes:"
  cat $OutDir/fimo_${ListLen}/fimo.tsv | grep "^$Motif" | cut -f3 | sort | uniq | grep -o -P "g.*?\.t\d" > $OutDir/fimo_${ListLen}/gene_containing_motifs.txt
  cat $OutDir/fimo_${ListLen}/gene_containing_motifs.txt | wc -l
done
```

## WC-1 Transcription factors

Motif were identified in the promoters of key groups of DEGs.

DEGs were first selected by:
* eye from heatmaps of DEGs under different conditions
* based on cluster analysis.


### Selected genes:

Attempt to identify common elements in manually selected groups of similarly regulated transcription factors.

```bash
screen -a
qlogin


ProjDir=/home/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics
cd $ProjDir
ProjDir=$(ls -d /home/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics)
Promoters=$(ls $ProjDir/analysis/promoters/V.dahliae/JR2/JR2_promoters.fa)
OutDir=analysis/promoters/motifs/V.dahliae/JR2/WC-1/transcription_factors/selected
mkdir -p $OutDir/meme

# Transcription factors each showing decresed expression in frq knock outs

GeneList="VDAG_JR2_Chr2g00510 VDAG_JR2_Chr2g08700 VDAG_JR2_Chr2g10990 VDAG_JR2_Chr5g01930 VDAG_JR2_Chr3g13030 VDAG_JR2_Chr2g10850 VDAG_JR2_Chr1g24030 VDAG_JR2_Chr3g06490 VDAG_JR2_Chr3g06100 VDAG_JR2_Chr2g03270 VDAG_JR2_Chr5g01890 VDAG_JR2_Chr8g07070 VDAG_JR2_Chr4g04620 VDAG_JR2_Chr8g07580 VDAG_JR2_Chr1g27370 VDAG_JR2_Chr1g03770 VDAG_JR2_Chr3g12060 VDAG_JR2_Chr2g00520 VDAG_JR2_Chr8g00200 VDAG_JR2_Chr1g26850 VDAG_JR2_Chr1g12350 VDAG_JR2_Chr2g10410 VDAG_JR2_Chr6g02190 VDAG_JR2_Chr7g07100 VDAG_JR2_Chr5g05720 VDAG_JR2_Chr2g01990 VDAG_JR2_Chr8g01630 VDAG_JR2_Chr5g00990 VDAG_JR2_Chr1g23950 VDAG_JR2_Chr2g07860 VDAG_JR2_Chr5g11380 VDAG_JR2_Chr2g07010 VDAG_JR2_Chr1g27300 VDAG_JR2_Chr5g05660 VDAG_JR2_Chr3g09510 VDAG_JR2_Chr2g00700 VDAG_JR2_Chr1g07150 VDAG_JR2_Chr3g08070 VDAG_JR2_Chr3g08730 VDAG_JR2_Chr2g04990
VDAG_JR2_Chr1g08930 VDAG_JR2_Chr1g08230 VDAG_JR2_Chr2g05130 VDAG_JR2_Chr3g04790 VDAG_JR2_Chr1g05210 VDAG_JR2_Chr7g10270 VDAG_JR2_Chr3g03860 VDAG_JR2_Chr8g11340 VDAG_JR2_Chr6g09990 VDAG_JR2_Chr1g09490 VDAG_JR2_Chr1g18820 VDAG_JR2_Chr2g11140 VDAG_JR2_Chr5g04520 VDAG_JR2_Chr6g07650 VDAG_JR2_Chr1g17310"
ListLen=$(echo $GeneList | grep -o 'g' | wc -l)
SubPromoters=$OutDir/promoters_${ListLen}.fa

for GeneID in $GeneList; do
cat $Promoters | sed -n "/^>${GeneID}_/,/^>/p" | grep -v "^$" | head -n -1
done > $SubPromoters

GeneList="VDAG_JR2_Chr6g09990 VDAG_JR2_Chr1g09490 VDAG_JR2_Chr1g18820 VDAG_JR2_Chr2g11140 VDAG_JR2_Chr5g04520 VDAG_JR2_Chr6g07650 VDAG_JR2_Chr1g17310"


ListLen=$(echo $GeneList | grep -o 'g' | wc -l)
SubPromoters=$OutDir/promoters_${ListLen}.fa

for GeneID in $GeneList; do
cat $Promoters | sed -n "/^>${GeneID}_/,/^>/p" | grep -v "^$" | head -n -1
done > $SubPromoters


meme $SubPromoters -dna -mod anr -nmotifs 5 -minw 6 -maxw 12 -revcomp -evt 1.0e+005 -oc $OutDir/meme_${ListLen}

cat $OutDir/meme_${ListLen}/meme.txt | grep -C2 'regular expression'

#---
# Running MAST
#---

mast $OutDir/meme_${ListLen}/meme.xml $SubPromoters -oc $OutDir/mast_${ListLen}
ls $OutDir/mast_${ListLen}/mast.txt

#---
# Running FIMO
#---
Thresh=0.00006 #CASSIS (trained for secmet clusters)
# Thresh=1e-4 #Default
fimo -thresh $Thresh -oc $OutDir/fimo_${ListLen} $OutDir/meme_${ListLen}/meme.html $OutDir/promoters_50.fa
for Motif in $(cat $OutDir/fimo_${ListLen}/fimo.tsv | grep 'TSS' | cut -f1 | sort | uniq); do
  echo $Motif
  echo "Promoters containing motif:"
  cat $OutDir/fimo_${ListLen}/fimo.tsv | grep "^$Motif" | cut -f3 | sort | uniq | wc -l
  echo "Total number found in these promoters:"
  cat $OutDir/fimo_${ListLen}/fimo.tsv | grep "^$Motif" | cut -f3 | sort | wc -l
  echo "this covers the following genes:"
  cat $OutDir/fimo_${ListLen}/fimo.tsv | grep "^$Motif" | cut -f3 | sort | uniq | grep -o -P "g.*?\.t\d" > $OutDir/fimo_${ListLen}/gene_containing_motifs.txt
  cat $OutDir/fimo_${ListLen}/gene_containing_motifs.txt | wc -l
done
```

```
ATTCDTTCATG
Promoters containing motif:
4
Total number found in these promoters:
5
this covers the following genes:
0
CAAGWGACA
Promoters containing motif:
11
Total number found in these promoters:
16
this covers the following genes:
0
CTKTCSACGAC
Promoters containing motif:
12
Total number found in these promoters:
17
this covers the following genes:
0
```


## Frq Transcription factors

Motif were identified in the promoters of key groups of DEGs.

DEGs were first selected by:
* eye from heatmaps of DEGs under different conditions
* based on cluster analysis.


### Selected genes:

Attempt to identify common elements in manually selected groups of similarly regulated transcription factors.

```bash
screen -a
qlogin


ProjDir=/home/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics
cd $ProjDir
# Promoters=$(ls out/g*.t1/PROMOTERS/all_promoter_sequences.fasta | head -n1)
ProjDir=$(ls -d /home/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics)
Promoters=$(ls $ProjDir/analysis/promoters/V.dahliae/JR2/JR2_promoters.fa)
OutDir=analysis/promoters/motifs/V.dahliae/JR2/transcription_factors/selected
mkdir -p $OutDir/meme

# Transcription factors each showing decresed expression in frq knock outs

GeneList="VDAG_JR2_Chr5g10370 VDAG_JR2_Chr8g11340 VDAG_JR2_Chr7g10270 VDAG_JR2_Chr2g11050 VDAG_JR2_Chr6g06230 VDAG_JR2_Chr4g09890 VDAG_JR2_Chr7g07100 VDAG_JR2_Chr5g05720 VDAG_JR2_Chr3g04200 VDAG_JR2_Chr2g07670 VDAG_JR2_Chr4g07710 VDAG_JR2_Chr1g05850 VDAG_JR2_Chr2g12610 VDAG_JR2_Chr6g08330 VDAG_JR2_Chr3g07150 VDAG_JR2_Chr3g08070 VDAG_JR2_Chr2g04990 VDAG_JR2_Chr2g11140 VDAG_JR2_Chr5g01890 VDAG_JR2_Chr1g23950 VDAG_JR2_Chr5g00990 VDAG_JR2_Chr7g02150 VDAG_JR2_Chr1g03770 VDAG_JR2_Chr2g13070 VDAG_JR2_Chr8g07070 VDAG_JR2_Chr2g00520 VDAG_JR2_Chr3g13030 VDAG_JR2_Chr2g00510 VDAG_JR2_Chr2g03990 VDAG_JR2_Chr5g10380 VDAG_JR2_Chr1g24030 VDAG_JR2_Chr8g00200 VDAG_JR2_Chr1g26850 VDAG_JR2_Chr2g10410 VDAG_JR2_Chr4g11420 VDAG_JR2_Chr2g10850 VDAG_JR2_Chr4g10010 VDAG_JR2_Chr8g01630 VDAG_JR2_Chr3g06100 VDAG_JR2_Chr2g03270 VDAG_JR2_Chr4g03340 VDAG_JR2_Chr8g08220 VDAG_JR2_Chr1g22480 VDAG_JR2_Chr1g03460 VDAG_JR2_Chr4g03300 VDAG_JR2_Chr5g05660 VDAG_JR2_Chr1g27300 VDAG_JR2_Chr1g10470 VDAG_JR2_Chr6g07650 VDAG_JR2_Chr5g05400"
ListLen=$(echo $GeneList | grep -o 'g' | wc -l)
SubPromoters=$OutDir/promoters_${ListLen}.fa

for GeneID in $GeneList; do
cat $Promoters | sed -n "/^>${GeneID}_/,/^>/p" | grep -v "^$" | head -n -1
done > $SubPromoters

# GeneList="VDAG_JR2_Chr5g10370 VDAG_JR2_Chr8g11340 VDAG_JR2_Chr7g10270 VDAG_JR2_Chr2g11050 VDAG_JR2_Chr6g06230"
GeneList="VDAG_JR2_Chr7g07100 VDAG_JR2_Chr5g05720 VDAG_JR2_Chr3g04200"
ListLen=$(echo $GeneList | grep -o 'g' | wc -l)
SubPromoters=$OutDir/promoters_${ListLen}.fa

for GeneID in $GeneList; do
cat $Promoters | sed -n "/^>${GeneID}_/,/^>/p" | grep -v "^$" | head -n -1
done > $SubPromoters


meme $SubPromoters -dna -mod anr -nmotifs 5 -minw 6 -maxw 12 -revcomp -evt 1.0e+005 -oc $OutDir/meme_${ListLen}

cat $OutDir/meme_${ListLen}/meme.txt | grep -C2 'regular expression'

#---
# Running MAST
#---

mast $OutDir/meme_${ListLen}/meme.xml $SubPromoters -oc $OutDir/mast_${ListLen}
ls $OutDir/mast_${ListLen}/mast.txt

#---
# Running FIMO
#---
Thresh=0.00006 #CASSIS (trained for secmet clusters)
# Thresh=1e-4 #Default
fimo -thresh $Thresh -oc $OutDir/fimo_${ListLen} $OutDir/meme_${ListLen}/meme.html $OutDir/promoters_50.fa
for Motif in $(cat $OutDir/fimo_${ListLen}/fimo.tsv | grep 'TSS' | cut -f1 | sort | uniq); do
  echo $Motif
  echo "Promoters containing motif:"
  cat $OutDir/fimo_${ListLen}/fimo.tsv | grep "^$Motif" | cut -f3 | sort | uniq | wc -l
  echo "Total number found in these promoters:"
  cat $OutDir/fimo_${ListLen}/fimo.tsv | grep "^$Motif" | cut -f3 | sort | wc -l
  echo "this covers the following genes:"
  cat $OutDir/fimo_${ListLen}/fimo.tsv | grep "^$Motif" | cut -f3 | sort | uniq | grep -o -P "g.*?\.t\d" > $OutDir/fimo_${ListLen}/gene_containing_motifs.txt
  cat $OutDir/fimo_${ListLen}/gene_containing_motifs.txt | wc -l
done
```

```
ATTCDTTCATG
Promoters containing motif:
4
Total number found in these promoters:
5
this covers the following genes:
0
CAAGWGACA
Promoters containing motif:
11
Total number found in these promoters:
16
this covers the following genes:
0
CTKTCSACGAC
Promoters containing motif:
12
Total number found in these promoters:
17
this covers the following genes:
0
```

# Identification using parameters optimised from N. crassa



## Iterations through different IG region lengths

```bash
# Assembly=$(ls ../clocks/public_genomes/JR2/Verticillium_dahliaejr2.GCA_000400815.2.dna.toplevel_parsed.fa)
Assembly=$(ls ../clocks/public_genomes/JR2/Verticillium_dahliaejr2.GCA_000400815.2.dna.toplevel.fa)
Gff=$(ls ../clocks/public_genomes/JR2/Verticillium_dahliaejr2.GCA_000400815.2.33.gff3)

OutDir=analysis/promoters/V.dahliae/JR2
mkdir -p $OutDir

for Dist in $(seq 1000 500 6000); do
  echo $Dist
  ProgDir=/home/armita/git_repos/emr_repos/scripts/verticillium_pathogenomics/analysis/promoter_analysis
  $ProgDir/extract_IG-regions_JR2.py --gff $Gff --fasta $Assembly --prefix $OutDir/JR2_IG_regions_${Dist} --distance $Dist
  # $ProgDir/extract_IG-regions_OR74.py --gff $Gff --fasta $Assembly --prefix $OutDir/OR74A_IG_regions_${Dist} --distance $Dist
  ls $OutDir/JR2_IG_regions_${Dist}*
done
```


```bash
cd /home/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics
OutDir=analysis/promoters/motifs/V.dahliae/JR2/WC-1/Smith2010
mkdir -p $OutDir

echo "VDAG_JR2_Chr1g01960
VDAG_JR2_Chr3g10380
VDAG_JR2_Chr1g17210
VDAG_JR2_Chr4g01610
VDAG_JR2_Chr3g07160
VDAG_JR2_Chr1g02000
VDAG_JR2_Chr1g14930
VDAG_JR2_Chr4g09240
VDAG_JR2_Chr6g10880
VDAG_JR2_Chr3g01670
VDAG_JR2_Chr6g00750
VDAG_JR2_Chr1g17240
VDAG_JR2_Chr3g10370
VDAG_JR2_Chr7g02530
VDAG_JR2_Chr1g17230" > $OutDir/gene_list_selected.txt

```


```bash
screen -a
qlogin

ProjDir=$(ls -d /home/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics)
cd $ProjDir

for Promoters in $(ls $ProjDir/analysis/promoters/V.dahliae/JR2/JR2_IG_regions_*.fa); do
  echo $Promoters
  PromLgth=$(echo $Promoters | rev | cut -f1 -d '_' | rev | sed 's/.fa//g')
  echo $PromLgth  
  OutDir=analysis/promoters/motifs/V.dahliae/JR2/WC-1/Smith2010/selected/${PromLgth}_intergenic_unmasked/meme
  mkdir -p $OutDir

  cat $Promoters | grep '>' | grep -v -f analysis/promoters/motifs/V.dahliae/JR2/WC-1/Smith2010/gene_list_selected.txt | tr -d '>' > $OutDir/gene_list_selected_control.txt
  cat $Promoters | grep '>' | grep 'kal-1' | tr -d '>'>> $OutDir/gene_list_selected_control.txt

  GeneList=$(cat analysis/promoters/motifs/V.dahliae/JR2/WC-1/Smith2010/gene_list_selected.txt)
  SubPromoters=$OutDir/promoters_OR74A_Smith2010_highpeaks.fa
  for GeneID in $GeneList; do
  cat $Promoters | sed -n "/[>|_]${GeneID}_/,/^>/p" | grep -v "^$" | head -n -1
  done > $SubPromoters

  GeneList=$(ls $OutDir/gene_list_selected_control.txt)
  ControlPromoters=$OutDir/promoters_OR74A_Smith2010_-ves.fa

  ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
  $ProgDir/extract_from_fasta.py --fasta $Promoters --headers $GeneList > $ControlPromoters
  ls $ControlPromoters
  # for GeneID in $GeneList; do
  # cat $Promoters | sed -n "/${GeneID}/,/^>/p" | grep -v "^$" | head -n -1
  # done > $ControlPromoters

  fasta-get-markov -m 2 $Promoters $OutDir/background_model.txt

  mkdir $OutDir/meme
  meme -objfun de -neg $ControlPromoters -bfile $OutDir/background_model.txt -dna -mod anr -nmotifs 5 -minw 4 -maxw 12 -revcomp -evt 0.05 -oc $OutDir/meme $SubPromoters
  # meme -objfun se -neg $ControlPromoters -bfile $OutDir/background_model.txt -dna -mod anr -nmotifs 5 -minw 4 -maxw 12 -revcomp -evt 0.05 -oc $OutDir/meme $SubPromoters
  # meme -neg -dna -mod anr -nmotifs 5 -minw 4 -maxw 12 -revcomp -evt 0.05 -oc $OutDir/meme $SubPromoters

  cat $OutDir/meme/meme.txt | grep -C2 'regular expression'

  #---
  # Running MAST
  #---

  mast $OutDir/meme/meme.xml $SubPromoters -oc $OutDir/mast
  # ls $OutDir/mast_${ListLen}/mast.txt


  #---
  # Running FIMO
  #---
  # Thresh=0.00006 #CASSIS (trained for secmet clusters)
  Thresh=1e-4 #Default
  fimo -thresh $Thresh -oc $OutDir/fimo $OutDir/meme/meme.html $Promoters
  for Motif in $(cat $OutDir/fimo/fimo.tsv | grep 'MEME-1' | cut -f1 | sort | uniq); do
    echo $Motif
    echo "Promoters containing motif:"
    cat $OutDir/fimo/fimo.tsv | grep "^$Motif" | cut -f3 | sort | uniq | wc -l
    echo "Total number found in these promoters:"
    cat $OutDir/fimo/fimo.tsv | grep "^$Motif" | cut -f3 | sort | wc -l
    echo "this covers the following genes:"
    cat $OutDir/fimo/fimo.tsv | grep "^$Motif" | cut -f3 | sort | uniq > $OutDir/fimo/gene_containing_motifs.txt
    cat $OutDir/fimo/gene_containing_motifs.txt | wc -l
  done
done 2>&1 | tee analysis/promoters/motifs/V.dahliae/JR2/WC-1/Smith2010/selected/intergenic_meme_unmasked.log

cat analysis/promoters/motifs/V.dahliae/JR2/WC-1/Smith2010/selected/intergenic_meme_unmasked.log | grep -e 'JR2_IG_regions_' -e 'Motif' | grep -B1 'Motif'| less -S
```

```
/home/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics/analysis/promoters/V.dahliae/JR2/JR2_IG_regions_5000.fa
	Motif AGGTASBTAGGK MEME-1 regular expression
  Promoters containing motif:
  3627
  Total number found in these promoters:
  10529
  this covers the following genes:
  3627
```



## Iterations through different promoter region lengths

```bash
# Assembly=$(ls ../clocks/public_genomes/JR2/Verticillium_dahliaejr2.GCA_000400815.2.dna.toplevel_parsed.fa)
Assembly=$(ls ../clocks/public_genomes/JR2/Verticillium_dahliaejr2.GCA_000400815.2.dna.toplevel.fa)
Gff=$(ls ../clocks/public_genomes/JR2/Verticillium_dahliaejr2.GCA_000400815.2.33.gff3)

OutDir=analysis/promoters/V.dahliae/JR2
mkdir -p $OutDir

for Dist in $(seq 1000 500 6000); do
  echo $Dist
  ProgDir=/home/armita/git_repos/emr_repos/scripts/verticillium_pathogenomics/analysis/promoter_analysis
  $ProgDir/extract_promoters_JR2.py --gff $Gff --fasta $Assembly --prefix $OutDir/JR2_promoters2_${Dist} --distance $Dist
  # $ProgDir/extract_IG-regions_OR74.py --gff $Gff --fasta $Assembly --prefix $OutDir/OR74A_IG_regions_${Dist} --distance $Dist
  ls $OutDir/JR2_promoters2_${Dist}*
done
```


```bash
cd /home/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics
OutDir=analysis/promoters/motifs/V.dahliae/JR2/WC-1/Smith2010
mkdir -p $OutDir

echo "VDAG_JR2_Chr1g01960
VDAG_JR2_Chr3g10380
VDAG_JR2_Chr1g17210
VDAG_JR2_Chr4g01610
VDAG_JR2_Chr3g07160
VDAG_JR2_Chr1g02000
VDAG_JR2_Chr1g14930
VDAG_JR2_Chr4g09240
VDAG_JR2_Chr6g10880
VDAG_JR2_Chr3g01670
VDAG_JR2_Chr6g00750
VDAG_JR2_Chr1g17240
VDAG_JR2_Chr3g10370
VDAG_JR2_Chr7g02530
VDAG_JR2_Chr1g17230" > $OutDir/gene_list_selected.txt

```


```bash
screen -a
qlogin

ProjDir=$(ls -d /home/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics)
cd $ProjDir

for Promoters in $(ls $ProjDir/analysis/promoters/V.dahliae/JR2/JR2_promoters2_*.fa); do
  echo $Promoters
  PromLgth=$(echo $Promoters | rev | cut -f1 -d '_' | rev | sed 's/.fa//g')
  echo $PromLgth  
  OutDir=analysis/promoters/motifs/V.dahliae/JR2/WC-1/Smith2010/selected/${PromLgth}_promoters_unmasked/meme
  mkdir -p $OutDir

  cat $Promoters | grep '>' | grep -v -f analysis/promoters/motifs/V.dahliae/JR2/WC-1/Smith2010/gene_list_selected.txt | tr -d '>' > $OutDir/gene_list_selected_control.txt
  # cat $Promoters | grep '>' | grep 'kal-1' | tr -d '>'>> $OutDir/gene_list_selected_control.txt

  # GeneList=$(cat analysis/promoters/motifs/V.dahliae/JR2/WC-1/Smith2010/gene_list_selected.txt)
  # SubPromoters=$OutDir/promoters_OR74A_Smith2010_highpeaks.fa
  # for GeneID in $GeneList; do
  # cat $Promoters | sed -n "/[>|_]${GeneID}_/,/^>/p" | grep -v "^$" | head -n -1
  # done > $SubPromoters


  cat $Promoters | grep '>' | grep -f analysis/promoters/motifs/V.dahliae/JR2/WC-1/Smith2010/gene_list_selected.txt | tr -d '>' > $OutDir/gene_list_selected_headers.txt

  GeneList=$(ls $OutDir/gene_list_selected_headers.txt)
  SubPromoters=$OutDir/promoters_OR74A_Smith2010_highpeaks.fa
  ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
  $ProgDir/extract_from_fasta.py --fasta $Promoters --headers $GeneList > $SubPromoters

  GeneList=$(ls $OutDir/gene_list_selected_control.txt)
  ControlPromoters=$OutDir/promoters_OR74A_Smith2010_-ves.fa

  ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
  $ProgDir/extract_from_fasta.py --fasta $Promoters --headers $GeneList > $ControlPromoters
  ls $ControlPromoters
  # for GeneID in $GeneList; do
  # cat $Promoters | sed -n "/${GeneID}/,/^>/p" | grep -v "^$" | head -n -1
  # done > $ControlPromoters

  fasta-get-markov -m 2 $Promoters $OutDir/background_model.txt

  mkdir $OutDir/meme
  meme -objfun de -neg $ControlPromoters -bfile $OutDir/background_model.txt -dna -mod anr -nmotifs 5 -minw 4 -maxw 12 -revcomp -evt 0.05 -oc $OutDir/meme $SubPromoters
  # meme -objfun se -neg $ControlPromoters -bfile $OutDir/background_model.txt -dna -mod anr -nmotifs 5 -minw 4 -maxw 12 -revcomp -evt 0.05 -oc $OutDir/meme $SubPromoters
  # meme -neg -dna -mod anr -nmotifs 5 -minw 4 -maxw 12 -revcomp -evt 0.05 -oc $OutDir/meme $SubPromoters

  cat $OutDir/meme/meme.txt | grep -C2 'regular expression'

  #---
  # Running MAST
  #---

  mast $OutDir/meme/meme.xml $SubPromoters -oc $OutDir/mast
  # ls $OutDir/mast_${ListLen}/mast.txt


  #---
  # Running FIMO
  #---
  # Thresh=0.00006 #CASSIS (trained for secmet clusters)
  Thresh=1e-4 #Default
  fimo -thresh $Thresh -oc $OutDir/fimo $OutDir/meme/meme.html $Promoters
  for Motif in $(cat $OutDir/fimo/fimo.tsv | grep 'MEME-' | cut -f1 | sort | uniq); do
    echo $Motif
    echo "Promoters containing motif:"
    cat $OutDir/fimo/fimo.tsv | grep "^$Motif" | cut -f3 | sort | uniq | wc -l
    echo "Total number found in these promoters:"
    cat $OutDir/fimo/fimo.tsv | grep "^$Motif" | cut -f3 | sort | wc -l
    echo "this covers the following genes:"
    cat $OutDir/fimo/fimo.tsv | grep "^$Motif" | cut -f3 | sort | uniq > $OutDir/fimo/gene_containing_motifs.txt
    cat $OutDir/fimo/gene_containing_motifs.txt | wc -l
  done
done 2>&1 | tee analysis/promoters/motifs/V.dahliae/JR2/WC-1/Smith2010/selected/promoters_meme_unmasked.log

cat analysis/promoters/motifs/V.dahliae/JR2/WC-1/Smith2010/selected/promoters_meme_unmasked.log | grep -e 'promoters_OR74A_Smith2010' -e 'Motif' | grep -B1 'Motif'| less -S
```

```
analysis/promoters/motifs/V.dahliae/JR2/WC-1/Smith2010/selected/2000_promoters_unmasked/meme/promoters_OR74A_Smith2010_-ves.fa
        Motif ACGT MEME-1 regular expression
```


This A[CG]GT motif was too short / common for mast and fimo to return significant results.


Enrichment of the motif to background sequences could be studied however using AME.

Note - Use of control promoters in this analysis showed that a hundreds / thousands of sequences
show more-significant enrichment of this sequence than the target sequences.

```bash
screen -a
qlogin

ProjDir=$(ls -d /home/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics)
cd $ProjDir

Promoters=$(ls analysis/promoters/V.dahliae/JR2/JR2_promoters2_2000.fa)
echo $Promoters
PromLgth=$(echo $Promoters | rev | cut -f1 -d '_' | rev | sed 's/.fa//g')
echo $PromLgth  
OutDir=analysis/promoters/motifs/V.dahliae/JR2/WC-1/Smith2010/selected/${PromLgth}_promoters_unmasked/meme
mkdir -p $OutDir

SubPromoters=$OutDir/promoters_OR74A_Smith2010_highpeaks.fa
ControlPromoters=$OutDir/promoters_OR74A_Smith2010_-ves.fa
Motif=$(ls $OutDir/meme/meme.txt)

ame --oc $OutDir/ame --kmer 2 --scoring avg --bfile $OutDir/background_model.txt $SubPromoters $Motif

cat $OutDir/ame/sequences.tsv  | cut -f3,6
#
# ame --oc $OutDir/ame --kmer 2 --scoring avg --bfile $OutDir/background_model.txt --control $ControlPromoters $SubPromoters $Motif
#
# cat $OutDir/ame/sequences.tsv  | cut -f3,6
#
# ame --oc $OutDir/ame_all --kmer 2 --scoring avg --bfile $OutDir/background_model.txt --control --shuffle-- $Promoters $Motif
```

Enrichment was detected in 8 of the sequences, with 2 being identified as false +ves:

```
seq_ID	class
VDAG_JR2_Chr1g01950_VDAG_JR2_Chr1g01960_merged_-2000	tp
VDAG_JR2_Chr1g17210_TSS_-2000_+50	tp
VDAG_JR2_Chr1g01990_VDAG_JR2_Chr1g02000_merged_-2000	tp
VDAG_JR2_Chr4g09240_TSS_-2000_+50	tp
VDAG_JR2_Chr6g10880_VDAG_JR2_Chr6g10890_merged_-2000	tp
VDAG_JR2_Chr3g01670_VDAG_JR2_Chr3g01680_merged_-2000	tp
VDAG_JR2_Chr1g17230_VDAG_JR2_Chr1g17240_merged_-2000	fp
VDAG_JR2_Chr1g17230_VDAG_JR2_Chr1g17240_merged_-2000	fp
```


MCAST was used to see if the motif is clustered


Note - mcast assumes that every time a binding site is observed it leads to binding.
(alpha = 1)

```bash
screen -a
qlogin

ProjDir=$(ls -d /home/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics)
cd $ProjDir

Promoters=$(ls analysis/promoters/V.dahliae/JR2/JR2_promoters2_2000.fa)
echo $Promoters
PromLgth=$(echo $Promoters | rev | cut -f1 -d '_' | rev | sed 's/.fa//g')
echo $PromLgth  
OutDir=analysis/promoters/motifs/V.dahliae/JR2/WC-1/Smith2010/selected/${PromLgth}_promoters_unmasked/meme
mkdir -p $OutDir

SubPromoters=$OutDir/promoters_OR74A_Smith2010_highpeaks.fa
ControlPromoters=$OutDir/promoters_OR74A_Smith2010_-ves.fa
Motif=$(ls $OutDir/meme/meme.txt)

mcast --alpha 0.31 --bfile $OutDir/background_model.txt --oc $OutDir/mcast $Motif $SubPromoters
```


# Iterations through promoter regions using RSAT

RSAT was found to give favourable results for N. crassa. Vd homologs as
identified above, were extracted from promoter regions and downloaded to my
local machine. These were uploaded to the RSAT webserver at:
http://rsat-tagc.univ-mrs.fr/oligo-analysis_form.cgi

running at oligomer lengths (kmers) 4-8, with a V. dahliae background model and
displaying on the webserver, rather than by email.

This was run for promoter lengths 1000 -> 6000bp at 500bp intervals. Results were
largely consistent accross promoter lengths as the sequences had many shared
promtoers, limiting the effect of length as these wouldnt change with sequence length.

Motifs consistently identified were:
GATCGATC	CAGGC	AGGGCTCC	ACAGCTCT	CGCTTGC	TGCGCA
using settings:
RSAT (words)	WCC	Smith2010 Vd promoter	15 sequences	6000bp length	unmasked	consensus 4-8	8bp Vd upstream no ORF		assembly of P<0.05 motifs.
```
$RSAT/perl-scripts/purge-sequence -i $RSAT/public_html/tmp/www-data/2019/09/11/tmp_sequence_2019-09-11.182035_zRJamW.fasta -format fasta -o $RSAT/public_html/tmp/www-data/2019/09/11/tmp_sequence_2019-09-11.182035_zRJamW.fasta.purged -seqtype dna;$RSAT/perl-scripts/oligo-analysis  -v 1 -sort -i $RSAT/public_html/tmp/www-data/2019/09/11/tmp_sequence_2019-09-11.182035_zRJamW.fasta.purged -format fasta  -lth occ_sig 0 -uth rank 50 -return occ,proba,rank -2str -noov -quick_if_possible  -seqtype dna -bg upstream-noorf -org Verticillium_dahliae_GCF_000150675.1_ASM15067v2 -pseudo 0.01 -l 4 -o $RSAT/public_html/tmp/www-data/2019/09/11/oligo-analysis_2019-09-11.182035_gHvuOd_4nt.tab; $RSAT/perl-scripts/oligo-analysis  -v 1 -sort -i $RSAT/public_html/tmp/www-data/2019/09/11/tmp_sequence_2019-09-11.182035_zRJamW.fasta.purged -format fasta  -lth occ_sig 0 -uth rank 50 -return occ,proba,rank -2str -noov -quick_if_possible  -seqtype dna -bg upstream-noorf -org Verticillium_dahliae_GCF_000150675.1_ASM15067v2 -pseudo 0.01 -l 5 -o $RSAT/public_html/tmp/www-data/2019/09/11/oligo-analysis_2019-09-11.182035_gHvuOd_5nt.tab; $RSAT/perl-scripts/oligo-analysis  -v 1 -sort -i $RSAT/public_html/tmp/www-data/2019/09/11/tmp_sequence_2019-09-11.182035_zRJamW.fasta.purged -format fasta  -lth occ_sig 0 -uth rank 50 -return occ,proba,rank -2str -noov -quick_if_possible  -seqtype dna -bg upstream-noorf -org Verticillium_dahliae_GCF_000150675.1_ASM15067v2 -pseudo 0.01 -l 6 -o $RSAT/public_html/tmp/www-data/2019/09/11/oligo-analysis_2019-09-11.182035_gHvuOd_6nt.tab; $RSAT/perl-scripts/oligo-analysis  -v 1 -sort -i $RSAT/public_html/tmp/www-data/2019/09/11/tmp_sequence_2019-09-11.182035_zRJamW.fasta.purged -format fasta  -lth occ_sig 0 -uth rank 50 -return occ,proba,rank -2str -noov -quick_if_possible  -seqtype dna -bg upstream-noorf -org Verticillium_dahliae_GCF_000150675.1_ASM15067v2 -pseudo 0.01 -l 7 -o $RSAT/public_html/tmp/www-data/2019/09/11/oligo-analysis_2019-09-11.182035_gHvuOd_7nt.tab; $RSAT/perl-scripts/oligo-analysis  -v 1 -sort -i $RSAT/public_html/tmp/www-data/2019/09/11/tmp_sequence_2019-09-11.182035_zRJamW.fasta.purged -format fasta  -lth occ_sig 0 -uth rank 50 -return occ,proba,rank -2str -noov -quick_if_possible  -seqtype dna -bg upstream-noorf -org Verticillium_dahliae_GCF_000150675.1_ASM15067v2 -pseudo 0.01 -l 8 -o $RSAT/public_html/tmp/www-data/2019/09/11/oligo-analysis_2019-09-11.182035_gHvuOd_8nt.tab; cat $RSAT/public_html/tmp/www-data/2019/09/11/oligo-analysis_2019-09-11.182035_gHvuOd_4nt.tab $RSAT/public_html/tmp/www-data/2019/09/11/oligo-analysis_2019-09-11.182035_gHvuOd_5nt.tab $RSAT/public_html/tmp/www-data/2019/09/11/oligo-analysis_2019-09-11.182035_gHvuOd_6nt.tab $RSAT/public_html/tmp/www-data/2019/09/11/oligo-analysis_2019-09-11.182035_gHvuOd_7nt.tab $RSAT/public_html/tmp/www-data/2019/09/11/oligo-analysis_2019-09-11.182035_gHvuOd_8nt.tab  -v 1 -sort -i $RSAT/public_html/tmp/www-data/2019/09/11/tmp_sequence_2019-09-11.182035_zRJamW.fasta.purged -format fasta  -lth occ_sig 0 -uth rank 50 -return occ,proba,rank -2str -noov -quick_if_possible  -seqtype dna -bg upstream-noorf -org Verticillium_dahliae_GCF_000150675.1_ASM15067v2 -pseudo 0.01
```

Results of the 1500bp and 6000bp promoter regions were downloaded to my local computer and
then copied onto the cluster:

```bash
scp -r /Users/armita/Downloads/aug_27th/Vd/RSAT/1500 cluster:/home/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics/analysis/promoters/motifs/V.dahliae/JR2/WC-1/Smith2010/selected/1500_promoters_unmasked/RSAT

scp -r /Users/armita/Downloads/aug_27th/Vd/RSAT/6000 cluster:/home/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics/analysis/promoters/motifs/V.dahliae/JR2/WC-1/Smith2010/selected/6000_promoters_unmasked/RSAT
```

From this, motifs were converted from transfac format to meme format and run on
original input sequences.


```bash
  cd /home/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics
  # for Dist in 1500 6000; do
  for Dist in 1500; do
    echo $Dist
    WorkDir=$(ls -d analysis/promoters/motifs/V.dahliae/JR2/WC-1/Smith2010/selected/${Dist}_promoters_unmasked/RSAT)
    Promoters=$(ls analysis/promoters/V.dahliae/JR2/JR2_promoters2_${Dist}.fa)
    SubPromoters=$(ls $WorkDir/../meme/promoters_OR74A_Smith2010_highpeaks.fa)
    Transfac=$(ls $WorkDir/*count_matrices.tf.txt)
    transfac2meme $Transfac > $WorkDir/motifs.meme

    # mkdir $WorkDir/mast
    # $WorkDir/fimo

    mast $WorkDir/motifs.meme $SubPromoters -oc $WorkDir/mast

    #---
    # Running FIMO
    #---
    # Thresh=0.00006 #CASSIS (trained for secmet clusters)
    Thresh=1e-4 #Default
    fimo -thresh $Thresh -oc $WorkDir/fimo $WorkDir/motifs.meme $Promoters
    for Motif in $(cat $WorkDir/fimo/fimo.tsv | grep 'assembly_' | cut -f1 | sort | uniq); do
    echo $Motif
    # echo "Promoters containing motif:"
    cat $WorkDir/fimo/fimo.tsv | grep "^$Motif" | cut -f3 | sort | uniq | wc -l
    # echo "Total number found in these promoters:"
    cat $WorkDir/fimo/fimo.tsv | grep "^$Motif" | cut -f3 | sort | wc -l
    # echo "this covers the following genes:"
    cat $WorkDir/fimo/fimo.tsv | grep "^$Motif" | cut -f3 | sort | uniq > $WorkDir/fimo/gene_containing_motifs.txt
    cat $WorkDir/fimo/gene_containing_motifs.txt | sed -r "s/_VDAG/\nVDAG/g" | wc -l
    echo
    done

    # mcast --bfile $WorkDir/../meme/background_model.txt --oc $WorkDir/mcast $WorkDir/motifs.meme $SubPromoters

    for Header in $(cat analysis/promoters/motifs/V.dahliae/JR2/WC-1/Smith2010/gene_list_selected.txt); do
    echo $Header
    cat $WorkDir/fimo/fimo.tsv | grep "$Header" | cut -f1 | sort | uniq -c | sed  "s/^ *//g"
    done
  done
```

Counts of FIMO results
Note that a palindromic sequence may appear at a site twice.
```
1500
assembly_1
1029
2288
1463

assembly_2
2891
3961
3998

assembly_3
1445
2132
1992

assembly_4
4083
7211
5669

assembly_5
2287
2851
3253

assembly_6
3261
4720
4504

assembly_7
1084
2394
1501


6000
assembly_1
1242
2914
1812

assembly_2
4396
9963
6272

assembly_3
3085
4922
4384

assembly_4
2730
3889
3971

assembly_5
3602
6338
5122

assembly_6
1341
3226
1917
```

Fimo gff annotations didnt import into geneious on my local computer properly
so I had to remove '+' and '-' characters from sequence names.

```bash
  cat /Users/armita/Downloads/aug_27th/Vd/selected/1500_promoters_unmasked/JR2_promoters2_1500.fa | tr -d '+-' > /Users/armita/Downloads/aug_27th/Vd/selected/1500_promoters_unmasked/JR2_promoters2_1500_edited.fa
  cat /Users/armita/Downloads/aug_27th/Vd/selected/1500_promoters_unmasked/RSAT/fimo/fimo.gff | sed 's/_[+-]/_/g' > /Users/armita/Downloads/aug_27th/Vd/selected/1500_promoters_unmasked/RSAT/fimo/fimo_edited.gff

  cat /Users/armita/Downloads/aug_27th/Vd/selected/6000_promoters_unmasked/JR2_promoters2_6000.fa | tr -d '+-' > /Users/armita/Downloads/aug_27th/Vd/selected/6000_promoters_unmasked/JR2_promoters2_6000_edited.fa
  cat /Users/armita/Downloads/aug_27th/Vd/selected/6000_promoters_unmasked/RSAT/fimo/fimo.gff | sed 's/_[+-]/_/g' > /Users/armita/Downloads/aug_27th/Vd/selected/6000_promoters_unmasked/RSAT/fimo/fimo_edited.gff
```

### Location of seven motifs from 1500bp promoters in 6000bp promoters:


```bash
  cd /home/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics

    WorkDir=$(ls -d analysis/promoters/motifs/V.dahliae/JR2/WC-1/Smith2010/selected/6000_promoters_unmasked/RSAT)
    Promoters=$(ls analysis/promoters/V.dahliae/JR2/JR2_promoters2_6000.fa)
    SubPromoters=$(ls $WorkDir/../meme/promoters_OR74A_Smith2010_highpeaks.fa)

    mast analysis/promoters/motifs/V.dahliae/JR2/WC-1/Smith2010/selected/1500_promoters_unmasked/RSAT/motifs.meme $SubPromoters -oc $WorkDir/mast_1500bp_motifs

    #---
    # Running FIMO
    #---
    # Thresh=0.00006 #CASSIS (trained for secmet clusters)
    Thresh=1e-4 #Default
    fimo -thresh $Thresh -oc $WorkDir/fimo_1500bp_motifs analysis/promoters/motifs/V.dahliae/JR2/WC-1/Smith2010/selected/1500_promoters_unmasked/RSAT/motifs.meme $Promoters
    for Motif in $(cat $WorkDir/fimo_1500bp_motifs/fimo.tsv | grep 'assembly_' | cut -f1 | sort | uniq); do
    echo $Motif
    # echo "Promoters containing motif:"
    cat $WorkDir/fimo_1500bp_motifs/fimo.tsv | grep "^$Motif" | cut -f3 | sort | uniq | wc -l
    # echo "Total number found in these promoters:"
    cat $WorkDir/fimo_1500bp_motifs/fimo.tsv | grep "^$Motif" | cut -f3 | sort | wc -l
    # echo "this covers the following genes:"
    cat $WorkDir/fimo_1500bp_motifs/fimo.tsv | grep "^$Motif" | cut -f3 | sort | uniq > $WorkDir/fimo_1500bp_motifs/gene_containing_motifs.txt
    cat $WorkDir/fimo_1500bp_motifs/gene_containing_motifs.txt | sed -r "s/_VDAG/\nVDAG/g" | wc -l
    echo
    done

    # mcast --bfile $WorkDir/../meme/background_model.txt --oc $WorkDir/mcast $WorkDir/motifs.meme $SubPromoters

    for Header in $(cat analysis/promoters/motifs/V.dahliae/JR2/WC-1/Smith2010/gene_list_selected.txt); do
    echo $Header
    cat $WorkDir/fimo_1500bp_motifs/fimo.tsv | grep "$Header" | cut -f1 | sort | uniq -c | sed  "s/^ *//g"
    done
```

```bash
cat /Users/armita/Downloads/aug_27th/Vd/selected/6000_promoters_unmasked/RSAT/fimo_1500bp_motifs/fimo.gff | sed 's/_[+-]/_/g' > /Users/armita/Downloads/aug_27th/Vd/selected/6000_promoters_unmasked/RSAT/fimo_1500bp_motifs/fimo_1500bp_edited.gff
```
