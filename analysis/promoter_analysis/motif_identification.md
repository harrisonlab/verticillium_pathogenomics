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

## Transcription factors

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
