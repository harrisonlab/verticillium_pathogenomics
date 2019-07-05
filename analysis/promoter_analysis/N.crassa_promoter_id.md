# Motif identification Nc

Commands used to identify motifs in the promoters of key groups of Nc WC-1 regulated genes

## Data download

The Nc genome OR74A (NC12) was downloaded:

```bash
cd /home/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics
OutDir=assembly/external_groups/N.crassa/OR74A
mkdir -p $OutDir
# data was downloaded to this location manually using my local machine
# from: https://www.ncbi.nlm.nih.gov/assembly/GCA_000182925.2
gunzip $OutDir/*.gz
```

Firstly, promoter sequences were extracted from the genome using
published gene models

```bash
Assembly=$(ls assembly/external_groups/N.crassa/OR74A/GCF_000182925.2_NC12_genomic.fna)
Gff=$(ls assembly/external_groups/N.crassa/OR74A/GCF_000182925.2_NC12_genomic.gff)

OutDir=analysis/promoters/N.crassa/OR74A
mkdir -p $OutDir
ProgDir=/home/armita/git_repos/emr_repos/scripts/verticillium_pathogenomics/analysis/promoter_analysis
for Dist in 500 1000 1500 2000 2500 3000; do
  echo $Dist
  $ProgDir/extract_promoters_OR74A.py --gff $Gff --fasta $Assembly --prefix $OutDir/OR74A_promoters_${Dist} --distance $Dist
done

ls $OutDir
```


## WC-1 characterised genes

Motif were identified in the promoters of key groups of DEGs.

DEGs were first selected by:
* eye from heatmaps of DEGs under different conditions
* based on cluster analysis.


A WC-1 motif has been identified in smith et al 2010.
Reidentification of this motif was attempted in the following genes:

```bash
  vvd, al-3, sub-1, csp-1, wc-1
```

These genes were identified in the OR74A gneome:

```
  vvd - gene8807
  al-3 - gene7756
  sub-1 - gene7260
  csp-1 - gene869
  wc-1 - gene9960
```

```bash
screen -a
qlogin


ProjDir=/home/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics
cd $ProjDir
ProjDir=$(ls -d /home/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics)
Promoters=$(ls $ProjDir/analysis/promoters/N.crassa/OR74A/OR74A_promoters.fa)
OutDir=analysis/promoters/motifs/N.crassa/OR74A/WC-1/selected_genes/selected
mkdir -p $OutDir/meme

GeneList="gene8807 gene7756 gene7260 gene869 gene9960"
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


## Smith et al selected genes


### Genes with highest ChIP-Seq peaks

#### unmasked Meme analysis

```bash
cd /home/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics
OutDir=analysis/promoters/motifs/N.crassa/OR74A/WC-1/Smith2010
mkdir -p $OutDir

echo "NCU02265 frq NCU03967 vvd NCU00582 cry NCU08699 bli-4 NCU03071 os-4/NCU03072 NCU02264/NCU02265 frqAS NCU10063 NCU04021 NCU05594 NCU02800/NCU02801 NCU07541 NCU00552 al-1 NCU03968 NCU11300/NCU06017 NCU00584/NCU00585 al-2" | grep -P -o "NCU\w+" > $OutDir/gene_list.txt

# printf $GeneList > $OutDir/gene_list.txt

Gff=$(ls assembly/external_groups/N.crassa/OR74A/GCF_000182925.2_NC12_genomic.gff)
cat $Gff | grep -P "\tgene\t" | grep -f $OutDir/gene_list.txt | cut -f9 | sed "s/;.*locus_tag=/\t/g" | cut -f1 -d ';' | sed 's/ID=//g' |  cut -f1 > $OutDir/gene_list_gene_id.txt
```


```bash
screen -a
qlogin

ProjDir=$(ls -d /home/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics)
cd $ProjDir

for Promoters in $(ls $ProjDir/analysis/promoters/N.crassa/OR74A/OR74A_promoters_*.fa); do
  echo $Promoters
  PromLgth=$(echo $Promoters | rev | cut -f1 -d '_' | rev | sed 's/.fa//g')
  echo $PromLgth  
  OutDir=analysis/promoters/motifs/N.crassa/OR74A/WC-1/Smith2010/selected/$PromLgth/meme2
  mkdir -p $OutDir

  GeneList=$(cat analysis/promoters/motifs/N.crassa/OR74A/WC-1/Smith2010/gene_list_gene_id.txt | sed "s/$/ /g" | tr -d '\n')
  SubPromoters=$OutDir/promoters_OR74A_Smith2010_highpeaks.fa
  for GeneID in $GeneList; do
  cat $Promoters | sed -n "/^>${GeneID}_/,/^>/p" | grep -v "^$" | head -n -1
  done > $SubPromoters

  mkdir $OutDir/meme
  meme $SubPromoters -dna -mod anr -nmotifs 5 -minw 6 -maxw 12 -revcomp -evt 0.05 -oc $OutDir/meme

  cat $OutDir/meme/meme.txt | grep -C2 'regular expression'

  #---
  # Running MAST
  #---

  mast $OutDir/meme/meme.xml $SubPromoters -oc $OutDir/mast
  # ls $OutDir/mast_${ListLen}/mast.txt
done > analysis/promoters/motifs/N.crassa/OR74A/WC-1/Smith2010/meme_unmasked2.log
# #---
# # Running FIMO
# #---
# # Thresh=0.00006 #CASSIS (trained for secmet clusters)
# Thresh=1e-4 #Default
# # fimo -thresh $Thresh -oc $OutDir/fimo_${ListLen} $OutDir/meme_${ListLen}/meme.html $OutDir/promoters_50.fa
# fimo -thresh $Thresh -oc $OutDir/fimo_${ListLen} $OutDir/meme_${ListLen}/meme.html $OutDir/promoters_${ListLen}.fa
# for Motif in $(cat $OutDir/fimo_${ListLen}/fimo.tsv | grep 'TSS' | cut -f1 | sort | uniq); do
#   echo $Motif
#   echo "Promoters containing motif:"
#   cat $OutDir/fimo_${ListLen}/fimo.tsv | grep "^$Motif" | cut -f3 | sort | uniq | wc -l
#   echo "Total number found in these promoters:"
#   cat $OutDir/fimo_${ListLen}/fimo.tsv | grep "^$Motif" | cut -f3 | sort | wc -l
#   echo "this covers the following genes:"
#   cat $OutDir/fimo_${ListLen}/fimo.tsv | grep "^$Motif" | cut -f3 | sort | uniq | grep -o -P "g.*?\.t\d" > $OutDir/fimo_${ListLen}/gene_containing_motifs.txt
#   cat $OutDir/fimo_${ListLen}/gene_containing_motifs.txt | wc -l
# done
```

<!--
```bash
screen -a
qlogin


ProjDir=/home/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics
cd $ProjDir
ProjDir=$(ls -d /home/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics)
Promoters=$(ls $ProjDir/analysis/promoters/N.crassa/OR74A/OR74A_promoters.fa)
OutDir=analysis/promoters/motifs/N.crassa/OR74A/WC-1/selected_genes/selected
mkdir -p $OutDir/dreme


OutDir=assembly/external_groups/N.crassa/OR74A_Smith2010
GeneList=$(cat $OutDir/gene_list_gene_id.txt | sed "s/$/ /g" | tr -d '\n')
ListLen=$(cat $OutDir/gene_list_gene_id.txt | wc -l)
SubPromoters=$OutDir/promoters_${ListLen}.fa

for GeneID in $GeneList; do
cat $Promoters | sed -n "/^>${GeneID}_/,/^>/p" | grep -v "^$" | head -n -1
done > $SubPromoters

# ControlPromoters=$OutDir/promoters_${ListLen}.fa
# for GeneID in $(cat $Promoters | grep '>' | tr -d '>' | cut -f1 -d '_' | grep -v -w -f $OutDir/gene_list_gene_id.txt); do
# cat $Promoters | sed -n "/^>${GeneID}_/,/^>/p" | grep -v "^$" | head -n -1
# done > $ControlPromoters

ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/promoters
$ProgDir/dreme_cut_promoters.py --fasta $SubPromoters > $OutDir/promoters_${ListLen}_fragmented_100bp.fa

dreme -verbosity 1 -oc $OutDir/dreme2 -dna -p $OutDir/promoters_${ListLen}_fragmented_100bp.fa -t 18000 -e 1 -m 5 -mink 6 -maxk 6
```

Dreme failed to identify the known motif


The analysis was repeated using 16 100bp long subsequences from WC-1 regulated genes that contained the known motif.

This file was made in Geneious and copied onto the cluster:

```bash
scp /Users/armita/Downloads/WC-1_100bp__ve_controls.fasta cluster:/home/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics/analysis/promoters/motifs/N.crassa/OR74A/WC-1/selected_genes/selected/.

scp /Users/armita/Downloads/WC-1_1000bp_pos_controls.fasta cluster:/home/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics/analysis/promoters/motifs/N.crassa/OR74A/WC-1/selected_genes/selected/.
```
 -->

```bash
screen -a
qlogin

ProjDir=$(ls -d /home/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics)
cd $ProjDir

for Promoters in $(ls $ProjDir/analysis/promoters/N.crassa/OR74A/OR74A_promoters_*.fa); do
  echo $Promoters
  PromLgth=$(echo $Promoters | rev | cut -f1 -d '_' | rev | sed 's/.fa//g')
  echo $PromLgth  
  OutDir=analysis/promoters/motifs/N.crassa/OR74A/WC-1/Smith2010/selected/$PromLgth/dreme_3
  mkdir -p $OutDir

  GeneList=$(cat analysis/promoters/motifs/N.crassa/OR74A/WC-1/Smith2010/gene_list_gene_id.txt | sed "s/$/ /g" | tr -d '\n')
  SubPromoters=$OutDir/promoters_OR74A_Smith2010_highpeaks.fa
  for GeneID in $GeneList; do
  cat $Promoters | sed -n "/^>${GeneID}_/,/^>/p" | grep -v "^$" | head -n -1
  done > $SubPromoters

  dreme -verbosity 1 -oc $OutDir/dreme -dna -p $SubPromoters -t 18000 -e 0.05 -mink 6 -maxk 12

  cat $OutDir/dreme/dreme.txt | grep 'MOTIF' | wc -l
  cat $OutDir/dreme/dreme.txt | grep 'GATCGA'
  #---
  # Running MAST
  #---
  mast $OutDir/dreme/dreme.xml $SubPromoters -oc $OutDir/mast
  # ls $OutDir/mast_${ListLen}/mast.txt
done > analysis/promoters/motifs/N.crassa/OR74A/WC-1/Smith2010/dreme_unmasked2.log
```

```bash
screen -a
qlogin

ProjDir=$(ls -d /home/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics)
cd $ProjDir

for Promoters in $(ls $ProjDir/analysis/promoters/N.crassa/OR74A/OR74A_promoters_*.fa); do
echo $Promoters
PromLgth=$(echo $Promoters | rev | cut -f1 -d '_' | rev | sed 's/.fa//g')
echo $PromLgth  
OutDir=analysis/promoters/motifs/N.crassa/OR74A/WC-1/Smith2010/selected/$PromLgth
mkdir -p $OutDir

GeneList=$(cat analysis/promoters/motifs/N.crassa/OR74A/WC-1/Smith2010/gene_list_gene_id.txt | sed "s/$/ /g" | tr -d '\n')
SubPromoters=$OutDir/promoters_OR74A_Smith2010_highpeaks.fa
for GeneID in $GeneList; do
cat $Promoters | sed -n "/^>${GeneID}_/,/^>/p" | grep -v "^$" | head -n -1
done > $SubPromoters

Better=0
MotifNum=0
while [[ $Better -lt 5 ]]; do
MotifNum=$(($MotifNum+1))
if [[ $MotifNum -gt 5 ]]; then
break
fi
echo $MotifNum
# Better=$(($Better+1))
mkdir -p $OutDir/glam/motif${MotifNum}
glam2 -O $OutDir/glam/motif${MotifNum} -a 6 -b 12 n $SubPromoters

cat $OutDir/glam/motif${MotifNum}/glam2.txt | sed -n '/Score/,/Ins/p;/Ins/q' | head -n -1
Score=$(cat $OutDir/glam/motif${MotifNum}/glam2.txt | sed -n '/Score/,/Ins/p;/Ins/q' | head -n1 | cut -f2 -d ' ')

mkdir -p $OutDir/glam/motif${MotifNum}/resample
printf "Best\t$Score\n" > $OutDir/glam/motif${MotifNum}/resample/scores.txt
for i in $(seq 1 100); do
shuffleseq -sequence $SubPromoters $OutDir/glam/motif${MotifNum}/resample/shuffled_${i}.fa
glam2 -O $OutDir/glam/motif${MotifNum}/resample/$i -a 6 -b 12 n $OutDir/glam/motif${MotifNum}/resample/shuffled_${i}.fa
ResampleScore=$(cat $OutDir/glam/motif${MotifNum}/resample/$i/glam2.txt | sed -n '/Score/,/Ins/p;/Ins/q' | head -n1 | cut -f2 -d ' ')
printf "$i\t$ResampleScore\n" >> $OutDir/glam/motif${MotifNum}/resample/scores.txt
done
Position=$(cat $OutDir/glam/motif${MotifNum}/resample/scores.txt | sort -k2 -n -r | grep -n 'Best' | cut -f1 -d ':')
Better=$(( $Position -1 ))
echo "${Better} motifs scored better than the GLAM2 motif${MotifNum}"
glam2mask -o $OutDir/glam/motif${MotifNum}/promoters_OR74A_Smith2010_highpeaks_masked${MotifNum}.fa $OutDir/glam/motif${MotifNum}/glam2.txt $SubPromoters
SubPromoters=$(ls $OutDir/glam/motif${MotifNum}/promoters_OR74A_Smith2010_highpeaks_masked${MotifNum}.fa)
done > $OutDir/glam2.log
done
```


### Smitch highestest peaks with repeatmasked genome


Promoter sequences were extracted from the masked genome using
published gene models

```bash
Assembly=$(ls assembly/external_groups/N.crassa/OR74A/GCF_000182925.2_NC12_genomic.fna)
MaskedAssembly=$(basename ${Assembly%.fna})_masked.fa
cat $Assembly | cut -f1 -d ' '| tr 'atgc' 'xxxx' > $MaskedAssembly
Gff=$(ls assembly/external_groups/N.crassa/OR74A/GCF_000182925.2_NC12_genomic.gff)

OutDir=analysis/promoters/N.crassa/OR74A
mkdir -p $OutDir
ProgDir=/home/armita/git_repos/emr_repos/scripts/verticillium_pathogenomics/analysis/promoter_analysis
for Dist in 500 1000 1500 2000 2500 3000; do
  echo $Dist
  $ProgDir/extract_promoters_OR74A.py --gff $Gff --fasta $MaskedAssembly --prefix $OutDir/OR74A_masked_promoters_${Dist} --distance $Dist
done

ls $OutDir
```

### Genes with highest ChIP-Seq peaks (masked)

#### Meme analysis (masked)

```bash
cd /home/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics
OutDir=analysis/promoters/motifs/N.crassa/OR74A/WC-1/Smith2010
mkdir -p $OutDir

echo "NCU02265 frq NCU03967 vvd NCU00582 cry NCU08699 bli-4 NCU03071 os-4/NCU03072 NCU02264/NCU02265 frqAS NCU10063 NCU04021 NCU05594 NCU02800/NCU02801 NCU07541 NCU00552 al-1 NCU03968 NCU11300/NCU06017 NCU00584/NCU00585 al-2" | grep -P -o "NCU\w+" > $OutDir/gene_list.txt

# printf $GeneList > $OutDir/gene_list.txt

Gff=$(ls assembly/external_groups/N.crassa/OR74A/GCF_000182925.2_NC12_genomic.gff)
cat $Gff | grep -P "\tgene\t" | grep -f $OutDir/gene_list.txt | cut -f9 | sed "s/;.*locus_tag=/\t/g" | cut -f1 -d ';' | sed 's/ID=//g' |  cut -f1 > $OutDir/gene_list_gene_id.txt
```


```bash
screen -a
qlogin

ProjDir=$(ls -d /home/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics)
cd $ProjDir

for Promoters in $(ls $ProjDir/analysis/promoters/N.crassa/OR74A/OR74A_masked_promoters_*.fa); do
  echo $Promoters
  PromLgth=$(echo $Promoters | rev | cut -f1 -d '_' | rev | sed 's/.fa//g')
  echo $PromLgth  
  OutDir=analysis/promoters/motifs/N.crassa/OR74A/WC-1/Smith2010/selected/${PromLgth}_masked/meme_2
  mkdir -p $OutDir

  GeneList=$(cat analysis/promoters/motifs/N.crassa/OR74A/WC-1/Smith2010/gene_list_gene_id.txt | sed "s/$/ /g" | tr -d '\n')
  SubPromoters=$OutDir/promoters_OR74A_Smith2010_highpeaks.fa
  for GeneID in $GeneList; do
  cat $Promoters | sed -n "/^>${GeneID}_/,/^>/p" | grep -v "^$" | head -n -1
  done > $SubPromoters

  mkdir $OutDir/meme
  meme $SubPromoters -dna -mod anr -nmotifs 5 -minw 6 -maxw 12 -revcomp -evt 0.05 -oc $OutDir/meme

  cat $OutDir/meme/meme.txt | grep -C2 'regular expression'

  #---
  # Running MAST
  #---

  mast $OutDir/meme/meme.xml $SubPromoters -oc $OutDir/mast
  # ls $OutDir/mast_${ListLen}/mast.txt
done > analysis/promoters/motifs/N.crassa/OR74A/WC-1/Smith2010/meme_masked2.log
```

### Dreme (masked)


```bash
screen -a
qlogin

ProjDir=$(ls -d /home/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics)
cd $ProjDir

for Promoters in $(ls $ProjDir/analysis/promoters/N.crassa/OR74A/OR74A_masked_promoters_*.fa); do
  echo $Promoters
  PromLgth=$(echo $Promoters | rev | cut -f1 -d '_' | rev | sed 's/.fa//g')
  echo $PromLgth  
  OutDir=analysis/promoters/motifs/N.crassa/OR74A/WC-1/Smith2010/selected/${PromLgth}_masked/dreme_2
  mkdir -p $OutDir

  GeneList=$(cat analysis/promoters/motifs/N.crassa/OR74A/WC-1/Smith2010/gene_list_gene_id.txt | sed "s/$/ /g" | tr -d '\n')
  SubPromoters=$OutDir/promoters_OR74A_Smith2010_highpeaks.fa
  for GeneID in $GeneList; do
  cat $Promoters | sed -n "/^>${GeneID}_/,/^>/p" | grep -v "^$" | head -n -1
  done > $SubPromoters

  # dreme -verbosity 1 -oc $OutDir/dreme -dna -p $SubPromoters -t 18000 -e 1.0e+005 -mink 6 -maxk 12
  dreme -verbosity 1 -oc $OutDir/dreme -dna -p $SubPromoters -t 18000 -e 0.05 -mink 6 -maxk 12


  cat $OutDir/dreme/dreme.txt | grep 'MOTIF' | wc -l
  cat $OutDir/dreme/dreme.txt | grep 'GATCGA'
  #---
  # Running MAST
  #---
  mast $OutDir/dreme/dreme.xml $SubPromoters -oc $OutDir/mast
  # ls $OutDir/mast_${ListLen}/mast.txt
done > analysis/promoters/motifs/N.crassa/OR74A/WC-1/Smith2010/dreme2.log
```

### Glam (masked)

```bash
screen -a
qlogin

ProjDir=$(ls -d /home/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics)
cd $ProjDir

for Promoters in $(ls $ProjDir/analysis/promoters/N.crassa/OR74A/OR74A_masked_promoters_*.fa); do
echo $Promoters
PromLgth=$(echo $Promoters | rev | cut -f1 -d '_' | rev | sed 's/.fa//g')
echo $PromLgth  
OutDir=analysis/promoters/motifs/N.crassa/OR74A/WC-1/Smith2010/selected/${PromLgth}_masked
mkdir -p $OutDir

GeneList=$(cat analysis/promoters/motifs/N.crassa/OR74A/WC-1/Smith2010/gene_list_gene_id.txt | sed "s/$/ /g" | tr -d '\n')
SubPromoters=$OutDir/promoters_OR74A_Smith2010_highpeaks.fa
for GeneID in $GeneList; do
cat $Promoters | sed -n "/^>${GeneID}_/,/^>/p" | grep -v "^$" | head -n -1
done > $SubPromoters

Better=0
MotifNum=0
while [[ $Better -lt 5 ]]; do
MotifNum=$(($MotifNum+1))
if [[ $MotifNum -gt 5 ]]; then
break
fi
echo $MotifNum
# Better=$(($Better+1))
mkdir -p $OutDir/glam/motif${MotifNum}
glam2 -O $OutDir/glam/motif${MotifNum} -a 6 -b 12 n $SubPromoters

cat $OutDir/glam/motif${MotifNum}/glam2.txt | sed -n '/Score/,/Ins/p;/Ins/q' | head -n -1
Score=$(cat $OutDir/glam/motif${MotifNum}/glam2.txt | sed -n '/Score/,/Ins/p;/Ins/q' | head -n1 | cut -f2 -d ' ')

mkdir -p $OutDir/glam/motif${MotifNum}/resample
printf "Best\t$Score\n" > $OutDir/glam/motif${MotifNum}/resample/scores.txt
for i in $(seq 1 100); do
shuffleseq -sequence $SubPromoters $OutDir/glam/motif${MotifNum}/resample/shuffled_${i}.fa
glam2 -O $OutDir/glam/motif${MotifNum}/resample/$i -a 6 -b 12 n $OutDir/glam/motif${MotifNum}/resample/shuffled_${i}.fa
ResampleScore=$(cat $OutDir/glam/motif${MotifNum}/resample/$i/glam2.txt | sed -n '/Score/,/Ins/p;/Ins/q' | head -n1 | cut -f2 -d ' ')
printf "$i\t$ResampleScore\n" >> $OutDir/glam/motif${MotifNum}/resample/scores.txt
done
Position=$(cat $OutDir/glam/motif${MotifNum}/resample/scores.txt | sort -k2 -n -r | grep -n 'Best' | cut -f1 -d ':')
Better=$(( $Position -1 ))
echo "${Better} motifs scored better than the GLAM2 motif${MotifNum}"
glam2mask -o $OutDir/glam/motif${MotifNum}/promoters_OR74A_Smith2010_highpeaks_masked${MotifNum}.fa $OutDir/glam/motif${MotifNum}/glam2.txt $SubPromoters
SubPromoters=$(ls $OutDir/glam/motif${MotifNum}/promoters_OR74A_Smith2010_highpeaks_masked${MotifNum}.fa)
done > $OutDir/glam2.log
done
```



## Defined +ve controls

```bash
screen -a
qlogin

ProjDir=/home/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics
cd $ProjDir
ProjDir=$(ls -d /home/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics)
OutDir=analysis/promoters/motifs/N.crassa/OR74A/WC-1/selected_genes/selected/dreme_controls
mkdir -p $OutDir

dreme -verbosity 1 -oc $OutDir -dna -p $OutDir/../WC-1_100bp__ve_controls.fasta -t 18000 -e 1 -m 5 -mink 6 -maxk 6

dreme -verbosity 1 -oc ${OutDir}_1000bp -dna -p $OutDir/../WC-1_1000bp_pos_controls.fasta -t 18000 -e 1 -m 5 -mink 6 -maxk 6

OutDir=analysis/promoters/motifs/N.crassa/OR74A/WC-1/selected_genes/selected/meme_controls
mkdir -p $OutDir
meme $OutDir/../WC-1_100bp__ve_controls.fasta -dna -mod anr -nmotifs 5 -minw 6 -maxw 12 -revcomp -evt 1.0e+005 -oc ${OutDir}_100bp
cat ${OutDir}_100bp/meme.txt | grep -C2 'regular expression'

meme $OutDir/../WC-1_1000bp_pos_controls.fasta -dna -mod anr -nmotifs 5 -minw 6 -maxw 12 -revcomp -evt 1.0e+005 -oc ${OutDir}_1000bp
cat ${OutDir}_1000bp/meme.txt | grep -C2 'regular expression'
```

These experiments show that the motif is identified in the 100bp +ve controls but not in the 1000bp +ve controls when using DREME and also when using MEME.

### False positive in literature?

Occurence of the GATCGA motif:

```bash
cat $Promoters | grep '>' | wc -l
# 9769
cat $Promoters | grep -n 'GATCGA' | wc -l
# 2323
cat $Promoters | grep -o -n 'GATCGA' | wc -l
# 2802
```

The high occurence of this motif may mean that this sequence is not specific to the WCC.

Emma identified that the reverse of this sequence (but not the reverse complement) is commonly present in the gapped motif that was identified in He and Liu 2005.


## He and Liu motif


### Glam analysis

```bash
cd /home/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics
OutDir=assembly/external_groups/N.crassa/OR74A_HeLiu2005
mkdir -p $OutDir
#
# echo "NCU02265 frq NCU03967 vvd NCU00582 cry NCU08699 bli-4 NCU03071 os-4/NCU03072 NCU02264/NCU02265 frqAS NCU10063 NCU04021 NCU05594 NCU02800/NCU02801 NCU07541 NCU00552 al-1 NCU03968 NCU11300/NCU06017 NCU00584/NCU00585 al-2" | grep -P -o "NCU\w+" > $OutDir/gene_list.txt
#
# # printf $GeneList > $OutDir/gene_list.txt
#
# Gff=$(ls assembly/external_groups/N.crassa/OR74A/GCF_000182925.2_NC12_genomic.gff)
# cat $Gff | grep -P "\tgene\t" | grep -f $OutDir/gene_list.txt | cut -f9 | sed "s/;.*locus_tag=/\t/g" | cut -f1 -d ';' | sed 's/ID=//g' |  cut -f1 > $OutDir/gene_list_gene_id.txt
```



<!-- Firstly, promoter sequences were extracted from the genome using
published gene models

```bash
Assembly=$(ls assembly/external_groups/N.crassa/OR74A/GCF_000182925.2_NC12_rna.fna)
Gff=$(ls assembly/external_groups/N.crassa/OR74A/GCF_000182925.2_NC12_genomic.gff)

OutDir=analysis/promoters/N.crassa/OR74A
mkdir -p $OutDir
ProgDir=/home/armita/git_repos/emr_repos/scripts/verticillium_pathogenomics/analysis/promoter_analysis
$ProgDir/extract_promoters_OR74A.py --gff $Gff --fasta $Assembly --prefix $OutDir/OR74A_promoters_rm --distance 2000

ls $OutDir
``` -->

```bash
screen -a
qlogin


ProjDir=/home/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics
cd $ProjDir
ProjDir=$(ls -d /home/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics)
Promoters=$(ls $ProjDir/analysis/promoters/N.crassa/OR74A/OR74A_promoters.fa)
OutDir=analysis/promoters/motifs/N.crassa/OR74A/WC-1/selected_genes/selected/glam
mkdir -p $OutDir

# GeneList=$(cat assembly/external_groups/N.crassa/OR74A_Smith2010/gene_list_gene_id.txt | sed "s/$/ /g" | tr -d '\n')
# ListLen=$(cat assembly/external_groups/N.crassa/OR74A_Smith2010/gene_list_gene_id.txt | wc -l)
# SubPromoters=$OutDir/promoters_${ListLen}.fa

#frq al-3 vvd
GeneList="gene10065 gene7756 gene8807"
ListLen=$(echo $GeneList | grep -o 'gene' | wc -l)
SubPromoters=$OutDir/promoters_${ListLen}.fa

# for GeneID in $GeneList; do
# cat $Promoters | sed -n "/^>${GeneID}_/,/^>/p" | grep -v "^$" | head -n -1
# done > $SubPromoters

# Promoters were masked
for GeneID in $GeneList; do
cat $Promoters | sed -n "/^>${GeneID}_/,/^>/p" | grep -v "^$" | head -n -1 | tr [:lower:] 'N'
done > $SubPromoters


# PurgedPromoters=$OutDir/promoters_${ListLen}_purged.fa
# purge -n -q -o $SubPromoters 1e-30 > $PurgedPromoters

glam2 -O $OutDir/glam_${ListLen} -a 10 -b 20 n $SubPromoters

glam2mask


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
# fimo -thresh $Thresh -oc $OutDir/fimo_${ListLen} $OutDir/meme_${ListLen}/meme.html $OutDir/promoters_50.fa
fimo -thresh $Thresh -oc $OutDir/fimo_${ListLen} $OutDir/meme_${ListLen}/meme.html $OutDir/promoters_${ListLen}.fa
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
