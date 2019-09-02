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

cat analysis/promoters/motifs/V.dahliae/JR2/WC-1/Smith2010/selected/intergenic_meme_unmasked.log | grep -e 'promoters_OR74A_Smith2010' -e 'Motif' | grep -B1 'Motif'| less -S
```

```
analysis/promoters/motifs/V.dahliae/JR2/WC-1/Smith2010/selected/4500_intergenic_unmasked/meme/promoters_OR74A_Smith2010_-ves.fa
        Motif TACTTAC MEME-1 regular expression

        Promoters containing motif:
        796
        Total number found in these promoters:
        936
        this covers the following genes:
        796
--
analysis/promoters/motifs/V.dahliae/JR2/WC-1/Smith2010/selected/6000_intergenic_unmasked/meme/promoters_OR74A_Smith2010_-ves.fa
        Motif GDGGWGGKRGG MEME-1 regular expression

        Promoters containing motif:
        5510
        Total number found in these promoters:
        26832
        this covers the following genes:
        5510
```
