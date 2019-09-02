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

# NCU07541 gene4665
# NCU10063 gene8039


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


# Smith manually curated promoters

N. crassa promoters were extracted for Smith 2010 genes.

Rather than using regions upstream of 5'UTRs, entire intergenic regions were selected. Furthermore some gene models were ignored, namely those where genes were nested in the 5'UTR of the Smith gene.


## RSAT

The extracted promoter set was first analysed by RSAT online with the following options:

* word analysis
* oligo anlaysis on kmer lengths 4-8bp
* 15 Smith et al promoter sequences
* comparison to N. crassa upstream regions (background model)

The 4th identified motif related to GATCGATC

```
AC  assembly_1
XX
ID  ttttattaataaaa
XX
DE  nnTTTTATTAATAAAAnn
P0       a     c     g     t
1        0     0     0     0
2        0     0     0     0
3        0     0     0 2.64577
4        0     0     0 4.86152
5        0     0     0 5.05831
6        0     0     0 6.18805
7    7.84257     0     0     0
8        0     0     0    10
9        0     0     0    10
10      10     0     0     0
11      10     0     0     0
12       0     0     0 7.84257
13   6.18805     0     0     0
14   5.05831     0     0     0
15   4.86152     0     0     0
16   2.64577     0     0     0
17       0     0     0     0
18       0     0     0     0
XX
CC  residue_type:
CC  program: transfac
CC  matrix.nb: 1
CC  accession: assembly_1
CC  AC: assembly_1
CC  id: assembly_1
CC  name: ttttattaataaaa
CC  version:
CC  name: ttttattaataaaa
CC  description: nnTTTTATTAATAAAAnn
CC  transfac_consensus:
CC  matrix.nb: 1
XX
//
AC  assembly_2
XX
ID  ttttattaaataataaaa
XX
DE  nnTTTTATTAAATAATAAAAnn
P0       a     c     g     t
1        0     0     0     0
2        0     0     0     0
3        0     0     0 2.64577
4        0     0     0 4.86152
5        0     0     0 5.05831
6        0     0     0 6.18805
7    7.84257     0     0     0
8        0     0     0    10
9        0     0     0    10
10      10     0     0     0
11      10     0     0     0
12   6.18805     0     0     0
13       0     0     0 6.18805
14   6.18805     0     0     0
15   5.51749     0     0     0
16       0     0     0 5.51749
17   5.51749     0     0     0
18   5.05831     0     0     0
19   2.64577     0     0     0
20   2.64577     0     0     0
21       0     0     0     0
22       0     0     0     0
XX
CC  residue_type:
CC  program: transfac
CC  matrix.nb: 2
CC  accession: assembly_2
CC  AC: assembly_2
CC  id: assembly_2
CC  name: ttttattaaataataaaa
CC  version:
CC  name: ttttattaaataataaaa
CC  description: nnTTTTATTAAATAATAAAAnn
CC  transfac_consensus:
CC  matrix.nb: 2
XX
//
AC  assembly_3
XX
ID  atactatagtattaataaaa
XX
DE  nnatACTATAGTATTAATAAAAnn
P0       a     c     g     t
1        0     0     0     0
2        0     0     0     0
3    0.0801749     0     0     0
4        0     0     0 0.0801749
5    1.17347     0     0     0
6        0 1.62536     0     0
7        0     0     0 4.26385
8    4.26385     0     0     0
9        0     0     0 4.26385
10   4.26385     0     0     0
11       0     0 1.62536     0
12       0     0     0 6.18805
13   7.84257     0     0     0
14       0     0     0    10
15       0     0     0    10
16      10     0     0     0
17      10     0     0     0
18       0     0     0 7.84257
19   6.18805     0     0     0
20   5.05831     0     0     0
21   4.86152     0     0     0
22   2.64577     0     0     0
23       0     0     0     0
24       0     0     0     0
XX
CC  residue_type:
CC  program: transfac
CC  matrix.nb: 3
CC  accession: assembly_3
CC  AC: assembly_3
CC  id: assembly_3
CC  name: atactatagtattaataaaa
CC  version:
CC  name: atactatagtattaataaaa
CC  description: nnatACTATAGTATTAATAAAAnn
CC  transfac_consensus:
CC  matrix.nb: 3
XX
//
AC  assembly_4
XX
ID  gatcgatc
XX
DE  nnGATCGATCnn
P0       a     c     g     t
1        0     0     0     0
2        0     0     0     0
3        0     0 3.73057     0
4       10     0     0     0
5        0     0     0    10
6        0    10     0     0
7        0     0    10     0
8       10     0     0     0
9        0     0     0    10
10       0 3.73057     0     0
11       0     0     0     0
12       0     0     0     0
XX
CC  residue_type:
CC  program: transfac
CC  matrix.nb: 4
CC  accession: assembly_4
CC  AC: assembly_4
CC  id: assembly_4
CC  name: gatcgatc
CC  version:
CC  name: gatcgatc
CC  description: nnGATCGATCnn
CC  transfac_consensus:
CC  matrix.nb: 4
XX
//
AC  assembly_5
XX
ID  cttaataaaa
XX
DE  nnCTTAATAAAAnn
P0       a     c     g     t
1        0     0     0     0
2        0     0     0     0
3        0 1.9898     0     0
4        0     0     0    10
5        0     0     0    10
6       10     0     0     0
7       10     0     0     0
8        0     0     0 7.84257
9    6.18805     0     0     0
10   5.05831     0     0     0
11   4.86152     0     0     0
12   2.64577     0     0     0
13       0     0     0     0
14       0     0     0     0
XX
CC  residue_type:
CC  program: transfac
CC  matrix.nb: 5
CC  accession: assembly_5
CC  AC: assembly_5
CC  id: assembly_5
CC  name: cttaataaaa
CC  version:
CC  name: cttaataaaa
CC  description: nnCTTAATAAAAnn
CC  transfac_consensus:
CC  matrix.nb: 5
XX
//
AC  assembly_6
XX
ID  gatcgatctc
XX
DE  nnGATCGATCTCnn
P0       a     c     g     t
1        0     0     0     0
2        0     0     0     0
3        0     0 3.73057     0
4       10     0     0     0
5        0     0     0    10
6        0    10     0     0
7        0     0    10     0
8       10     0     0     0
9        0     0     0    10
10       0 5.90674     0     0
11       0     0     0 5.90674
12       0 5.90674     0     0
13       0     0     0     0
14       0     0     0     0
XX
CC  residue_type:
CC  program: transfac
CC  matrix.nb: 6
CC  accession: assembly_6
CC  AC: assembly_6
CC  id: assembly_6
CC  name: gatcgatctc
CC  version:
CC  name: gatcgatctc
CC  description: nnGATCGATCTCnn
CC  transfac_consensus:
CC  matrix.nb: 6
XX
//
AC  assembly_7
XX
ID  ttttattaattaaataaaa
XX
DE  nnTTTTATTAATTAAATAAAAnn
P0       a     c     g     t
1        0     0     0     0
2        0     0     0     0
3        0     0     0 2.64577
4        0     0     0 4.86152
5        0     0     0 5.05831
6        0     0     0 6.18805
7    7.84257     0     0     0
8        0     0     0    10
9        0     0     0    10
10      10     0     0     0
11      10     0     0     0
12       0     0     0 7.84257
13       0     0     0 5.05831
14   5.05831     0     0     0
15   6.18805     0     0     0
16   6.18805     0     0     0
17       0     0     0 6.18805
18   6.18805     0     0     0
19   2.64577     0     0     0
20   2.64577     0     0     0
21   2.64577     0     0     0
22       0     0     0     0
23       0     0     0     0
XX
CC  residue_type:
CC  program: transfac
CC  matrix.nb: 7
CC  accession: assembly_7
CC  AC: assembly_7
CC  id: assembly_7
CC  name: ttttattaattaaataaaa
CC  version:
CC  name: ttttattaattaaataaaa
CC  description: nnTTTTATTAATTAAATAAAAnn
CC  transfac_consensus:
CC  matrix.nb: 7
XX
//
AC  assembly_8
XX
ID  ttttattaatactataaaa
XX
DE  nnTTTTATTAATACTATAAAAnn
P0       a     c     g     t
1        0     0     0     0
2        0     0     0     0
3        0     0     0 2.64577
4        0     0     0 4.86152
5        0     0     0 5.05831
6        0     0     0 6.18805
7    7.84257     0     0     0
8        0     0     0    10
9        0     0     0    10
10      10     0     0     0
11      10     0     0     0
12       0     0     0 7.84257
13   6.18805     0     0     0
14       0 1.62536     0     0
15       0     0     0 4.26385
16   4.26385     0     0     0
17       0     0     0 4.86152
18   4.86152     0     0     0
19   4.86152     0     0     0
20   4.86152     0     0     0
21   2.64577     0     0     0
22       0     0     0     0
23       0     0     0     0
XX
CC  residue_type:
CC  program: transfac
CC  matrix.nb: 8
CC  accession: assembly_8
CC  AC: assembly_8
CC  id: assembly_8
CC  name: ttttattaatactataaaa
CC  version:
CC  name: ttttattaatactataaaa
CC  description: nnTTTTATTAATACTATAAAAnn
CC  transfac_consensus:
CC  matrix.nb: 8
XX
//
AC  assembly_9
XX
ID  ttttattaatatagtat
XX
DE  nnTTTTATTAATATAGTatnn
P0       a     c     g     t
1        0     0     0     0
2        0     0     0     0
3        0     0     0 2.64577
4        0     0     0 4.86152
5        0     0     0 5.05831
6        0     0     0 6.18805
7    7.84257     0     0     0
8        0     0     0    10
9        0     0     0    10
10      10     0     0     0
11      10     0     0     0
12       0     0     0 7.84257
13   6.18805     0     0     0
14       0     0     0 4.26385
15   4.26385     0     0     0
16       0     0 1.62536     0
17       0     0     0 1.17347
18   0.0801749     0     0     0
19       0     0     0 0.0801749
20       0     0     0     0
21       0     0     0     0
XX
CC  residue_type:
CC  program: transfac
CC  matrix.nb: 9
CC  accession: assembly_9
CC  AC: assembly_9
CC  id: assembly_9
CC  name: ttttattaatatagtat
CC  version:
CC  name: ttttattaatatagtat
CC  description: nnTTTTATTAATATAGTatnn
CC  transfac_consensus:
CC  matrix.nb: 9
XX
//
AC  assembly_10
XX
ID  ttttattaatagtatag
XX
DE  nnTTTTATTAATAGTATAGnn
P0       a     c     g     t
1        0     0     0     0
2        0     0     0     0
3        0     0     0 2.64577
4        0     0     0 4.86152
5        0     0     0 5.05831
6        0     0     0 6.18805
7    7.84257     0     0     0
8        0     0     0    10
9        0     0     0    10
10      10     0     0     0
11      10     0     0     0
12       0     0     0 7.84257
13   6.18805     0     0     0
14       0     0 1.17347     0
15       0     0     0 4.26385
16   4.26385     0     0     0
17       0     0     0 4.26385
18   4.26385     0     0     0
19       0     0 1.62536     0
20       0     0     0     0
21       0     0     0     0
XX
CC  residue_type:
CC  program: transfac
CC  matrix.nb: 10
CC  accession: assembly_10
CC  AC: assembly_10
CC  id: assembly_10
CC  name: ttttattaatagtatag
CC  version:
CC  name: ttttattaatagtatag
CC  description: nnTTTTATTAATAGTATAGnn
CC  transfac_consensus:
CC  matrix.nb: 10
XX
//
```

The 4th motif was extracted from this output and saved as a text file, which could be converted into meme format using transfec2motif.pl


(From local machine)
```bash
./scripts/transfac2meme /Users/armita/Documents/GATCGATC_motif_transfec.txt > /Users/armita/Documents/GATCGATC_motif_transfec.meme
```

This was uploaded onto the cluster, along with the original fasta file where Fimo was run to assess occurence in the original sequences:

```bash
scp /Users/armita/Downloads/Smith\ promoters\ untrimmed.fasta cluster:/home/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics/analysis/promoters/N.crassa/OR74A/RSAT/smith_promoters_untrimmed.fasta
```

```bash
cat /home/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics/analysis/promoters/N.crassa/OR74A/RSAT/smith_promoters_untrimmed.fasta | sed 's/ /_/g' > /home/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics/analysis/promoters/N.crassa/OR74A/RSAT/smith_promoters_untrimmed_parsed.fasta

fimo -o /home/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics/analysis/promoters/N.crassa/OR74A/RSAT/fimo_out /home/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics/analysis/promoters/N.crassa/OR74A/RSAT/GATCGATC_motif_transfec.meme /home/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics/analysis/promoters/N.crassa/OR74A/RSAT/smith_promoters_untrimmed_parsed.fasta

cat /home/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics/analysis/promoters/N.crassa/OR74A/RSAT/fimo_out/fimo.tsv | tail -n+2 | grep -v '#' | sed -r "s/\s+/ /g" | cut -f2 -d ' ' | sort | uniq -c | sort -nr
```

```
12 NCU03968_promoter
 8 vvd_NCU03967_promoter_(reversed)
 8 OS-4_NCU03071_and_NCU03072_promoter
 6 cry-dash_NCU00582_promoter
 4 NCU10063_promoter_(reversed)
 4 frq_NCU02264_promoter
 2 NCU04021_promoter
 2 NCU02800_and_NCU02801_promoter
 2 bli-4_NCU08699_promoter
```


## Running using MEME

Due to the success of RSAT, the same results were attempted to be identified using MEME.

```bash
screen -a
qlogin

ProjDir=$(ls -d /home/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics)
cd $ProjDir

for Promoters in $(ls $ProjDir/analysis/promoters/N.crassa/OR74A/RSAT/smith_promoters_untrimmed_parsed.fasta); do
  echo $Promoters
  PromLgth="intergenic"
  echo $PromLgth  
  OutDir=analysis/promoters/motifs/N.crassa/OR74A/WC-1/Smith2010/selected/${PromLgth}_unmasked/meme
  mkdir -p $OutDir

  mkdir $OutDir/meme
  meme $Promoters -dna -mod anr -nmotifs 20 -minw 4 -maxw 12 -revcomp -evt 0.05 -oc $OutDir/meme

  cat $OutDir/meme/meme.txt | grep -C2 'regular expression'

  #---
  # Running MAST
  #---

  mast $OutDir/meme/meme.xml $SubPromoters -oc $OutDir/mast
done
```


# Running using total Intergenic regions


```bash
Assembly=$(ls assembly/external_groups/N.crassa/OR74A/GCF_000182925.2_NC12_genomic.fna)
Gff=$(ls assembly/external_groups/N.crassa/OR74A/GCF_000182925.2_NC12_genomic.gff)

OutDir=analysis/promoters/N.crassa/OR74A
mkdir -p $OutDir
Dist=5000
# ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/promoters
# $ProgDir/extract_IG-regions.py --gff $Gff --fasta $Assembly --prefix $OutDir/OR74A_IG_regions
ProgDir=/home/armita/git_repos/emr_repos/scripts/verticillium_pathogenomics/analysis/promoter_analysis
$ProgDir/extract_IG-regions_OR74.py --gff $Gff --fasta $Assembly --prefix $OutDir/OR74A_IG_regions_${Dist} --distance $Dist

ls $OutDir/OR74A_IG_regions_${Dist}*
```


```bash
cd /home/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics
OutDir=analysis/promoters/motifs/N.crassa/OR74A/WC-1/Smith2010
mkdir -p $OutDir

echo "frq vvd cry bli-4 os-4 NCU02264 NCU10063 NCU04021 NCU05594 NCU02800 NCU07541 al-1 NCU03968 NCU11300 al-2" | sed 's/ /\n/g' > $OutDir/gene_list_selected.txt

Promoters=$(ls analysis/promoters/N.crassa/OR74A/OR74A_IG_regions_5000.fa)
cat $Promoters | grep '>' | grep -v -f $OutDir/gene_list_selected.txt | tr -d '>' > $OutDir/gene_list_selected_control.txt
cat $Promoters | grep '>' | grep 'kal-1' | tr -d '>'>> $OutDir/gene_list_selected_control.txt

```


```bash
screen -a
qlogin

ProjDir=$(ls -d /home/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics)
cd $ProjDir

for Promoters in $(ls $ProjDir/analysis/promoters/N.crassa/OR74A/OR74A_IG_regions_5000.fa); do
  echo $Promoters
  PromLgth=$(echo $Promoters | rev | cut -f1 -d '_' | rev | sed 's/.fa//g')
  echo $PromLgth  
  OutDir=analysis/promoters/motifs/N.crassa/OR74A/WC-1/Smith2010/selected/${PromLgth}_intergenic_unmasked/meme
  mkdir -p $OutDir

  GeneList=$(cat analysis/promoters/motifs/N.crassa/OR74A/WC-1/Smith2010/gene_list_selected.txt)
  SubPromoters=$OutDir/promoters_OR74A_Smith2010_highpeaks.fa
  for GeneID in $GeneList; do
  cat $Promoters | sed -n "/[>|_]${GeneID}_/,/^>/p" | grep -v "^$" | head -n -1
  done > $SubPromoters

  GeneList=$(ls analysis/promoters/motifs/N.crassa/OR74A/WC-1/Smith2010/gene_list_selected_control.txt)
  ControlPromoters=$OutDir/promoters_OR74A_Smith2010_-ves.fa

  ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
  $ProgDir/extract_from_fasta.py --fasta $Promoters --headers $GeneList > $ControlPromoters
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
done > analysis/promoters/motifs/N.crassa/OR74A/WC-1/Smith2010/meme_masked.log
```

# Smith et al 2010 TF genes


```bash
cd /home/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics
OutDir=analysis/promoters/motifs/N.crassa/OR74A/WC-1/Smith2010
mkdir -p $OutDir

echo "sah-1
wc-1
NCU04295
NCU09615
nit-2
NCU05964
NCU09829
NCU01243
NCU01871
hsf-2
csp-1
cre-1
adv-1
bek-1
NCU07705
sub-1
NCU06534
sah-2
sre
NCU07846
NCU00275
vad-2
NCU05994
NCU03273
NCU03184
NCU08000
NCU08159
NCU06095" | sed 's/ /\n/g' > $OutDir/gene_list_TF_selected.txt

Promoters=$(ls analysis/promoters/N.crassa/OR74A/OR74A_IG_regions_5000.fa)
cat $Promoters | grep '>' | grep -v -f $OutDir/gene_list_TF_selected.txt | tr -d '>' > $OutDir/gene_list_TF_selected_control.txt
# cat $Promoters | grep '>' | grep 'kal-1' | tr -d '>'>> $OutDir/gene_list_TF_selected_control.txt
```


```bash
screen -a
qlogin

ProjDir=$(ls -d /home/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics)
cd $ProjDir

for Promoters in $(ls $ProjDir/analysis/promoters/N.crassa/OR74A/OR74A_IG_regions_5000.fa); do
  echo $Promoters
  PromLgth=$(echo $Promoters | rev | cut -f1 -d '_' | rev | sed 's/.fa//g')
  echo $PromLgth  
  OutDir=analysis/promoters/motifs/N.crassa/OR74A/WC-1/Smith2010/selected_TF/${PromLgth}_intergenic_unmasked/meme
  mkdir -p $OutDir

  GeneList=$(cat analysis/promoters/motifs/N.crassa/OR74A/WC-1/Smith2010/gene_list_TF_selected.txt)
  SubPromoters=$OutDir/promoters_OR74A_Smith2010_highpeaks.fa
  for GeneID in $GeneList; do
  cat $Promoters | sed -n "/[>|_]${GeneID}_/,/^>/p" | grep -v "^$" | head -n -1
  done > $SubPromoters

  GeneList=$(ls analysis/promoters/motifs/N.crassa/OR74A/WC-1/Smith2010/gene_list_TF_selected_control.txt)
  ControlPromoters=$OutDir/promoters_OR74A_Smith2010_-ves.fa

  ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
  $ProgDir/extract_from_fasta.py --fasta $Promoters --headers $GeneList > $ControlPromoters
  # for GeneID in $GeneList; do
  # cat $Promoters | sed -n "/${GeneID}/,/^>/p" | grep -v "^$" | head -n -1
  # done > $ControlPromoters

  fasta-get-markov -m 2 $Promoters $OutDir/background_model.txt

  mkdir $OutDir/meme
  # meme -objfun de -neg $ControlPromoters -bfile $OutDir/background_model.txt -dna -mod anr -nmotifs 5 -minw 4 -maxw 12 -revcomp -evt 0.05 -oc $OutDir/meme $SubPromoters
  meme -objfun se -neg $ControlPromoters -bfile $OutDir/background_model.txt -dna -mod anr -nmotifs 5 -minw 4 -maxw 12 -revcomp -evt 0.05 -oc $OutDir/meme2 $SubPromoters
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
done > analysis/promoters/motifs/N.crassa/OR74A/WC-1/Smith2010/meme_masked.log
```


# Smith et al 2010 TF genes subselection


```bash
cd /home/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics
OutDir=analysis/promoters/motifs/N.crassa/OR74A/WC-1/Smith2010
mkdir -p $OutDir

echo "sah-1
wc-1
NCU04295
NCU09615
nit-2
NCU05964
NCU09829
NCU01243
NCU01871" | sed 's/ /\n/g' > $OutDir/gene_list_TF2_selected.txt

Promoters=$(ls analysis/promoters/N.crassa/OR74A/OR74A_IG_regions_5000.fa)
cat $Promoters | grep '>' | grep -v -f $OutDir/gene_list_TF2_selected.txt | tr -d '>' > $OutDir/gene_list_TF2_selected_control.txt

```


```bash
screen -a
qlogin

ProjDir=$(ls -d /home/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics)
cd $ProjDir

for Promoters in $(ls $ProjDir/analysis/promoters/N.crassa/OR74A/OR74A_IG_regions_5000.fa); do
  echo $Promoters
  PromLgth=$(echo $Promoters | rev | cut -f1 -d '_' | rev | sed 's/.fa//g')
  echo $PromLgth  
  OutDir=analysis/promoters/motifs/N.crassa/OR74A/WC-1/Smith2010/selected_TF2/${PromLgth}_intergenic_unmasked/meme
  mkdir -p $OutDir

  GeneList=$(cat analysis/promoters/motifs/N.crassa/OR74A/WC-1/Smith2010/gene_list_TF2_selected.txt)
  SubPromoters=$OutDir/promoters_OR74A_Smith2010_highpeaks.fa
  for GeneID in $GeneList; do
  cat $Promoters | sed -n "/[>|_]${GeneID}_/,/^>/p" | grep -v "^$" | head -n -1
  done > $SubPromoters

  GeneList=$(ls analysis/promoters/motifs/N.crassa/OR74A/WC-1/Smith2010/gene_list_TF2_selected_control.txt)
  ControlPromoters=$OutDir/promoters_OR74A_Smith2010_-ves.fa

  ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
  $ProgDir/extract_from_fasta.py --fasta $Promoters --headers $GeneList > $ControlPromoters
  # for GeneID in $GeneList; do
  # cat $Promoters | sed -n "/${GeneID}/,/^>/p" | grep -v "^$" | head -n -1
  # done > $ControlPromoters

  fasta-get-markov -m 2 $ControlPromoters $OutDir/background_model.txt

  mkdir $OutDir/meme
  meme -objfun de -neg $ControlPromoters -bfile $OutDir/background_model.txt -dna -mod anr -nmotifs 5 -minw 4 -maxw 12 -revcomp -evt 0.05 -oc $OutDir/meme $SubPromoters
  # meme -objfun se -neg $ControlPromoters -bfile $OutDir/background_model.txt -dna -mod anr -nmotifs 5 -minw 4 -maxw 12 -revcomp -evt 0.05 -oc $OutDir/meme2 $SubPromoters
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
done
```

This identified 2838 genes, contained 4151 occurences of the motif


## Searches for GATCGA motif in TF genes

```bash
Motif=$(ls analysis/promoters/motifs/N.crassa/OR74A/WC-1/Smith2010/selected/5000_intergenic_unmasked/meme/meme/meme.html)
Promoters=$(ls analysis/promoters/motifs/N.crassa/OR74A/WC-1/Smith2010/selected_TF/5000_intergenic_unmasked/meme/promoters_OR74A_Smith2010_highpeaks.fa)

OutDir=analysis/promoters/motifs/N.crassa/OR74A/WC-1/Smith2010/selected/5000_intergenic_unmasked/meme
Thresh=1e-4 #Default
fimo -thresh $Thresh -oc $OutDir/fimoTF $Motif $Promoters
for Motif in $(cat $OutDir/fimo/fimo.tsv | grep 'MEME-1' | cut -f1 | sort | uniq); do
  echo $Motif
  echo "Promoters containing motif:"
  cat $OutDir/fimoTF/fimo.tsv | grep "^$Motif" | cut -f3 | sort | uniq | wc -l
  echo "Total number found in these promoters:"
  cat $OutDir/fimoTF/fimo.tsv | grep "^$Motif" | cut -f3 | sort | wc -l
  echo "this covers the following genes:"
  cat $OutDir/fimoTF/fimo.tsv | grep "^$Motif" | cut -f3 | sort | uniq > $OutDir/fimoTF/gene_containing_motifs.txt
  cat $OutDir/fimoTF/gene_containing_motifs.txt | wc -l
done
```


18 occurences were found in 8 genes. This represented 8 of the top 9 TF genes, as ranked by z-score in the smith paper.
```
mig-12_NCU09829_IG_region-5000
NCU01243_IG_region-5000
NCU01871_NCU01872_IG_region-5000
NCU05964_NCU17127_IG_region-5000
NCU09615_IG_region-5000
nit-2_NCU09066_IG_region-5000
sah-1_IG_region-5000
wc-1_IG_region-5000
```


## Iterations through different IG region lengths



```bash
Assembly=$(ls assembly/external_groups/N.crassa/OR74A/GCF_000182925.2_NC12_genomic.fna)
Gff=$(ls assembly/external_groups/N.crassa/OR74A/GCF_000182925.2_NC12_genomic.gff)

OutDir=analysis/promoters/N.crassa/OR74A
mkdir -p $OutDir
for Dist in $(seq 1000 500 6000); do
  echo $Dist
  ProgDir=/home/armita/git_repos/emr_repos/scripts/verticillium_pathogenomics/analysis/promoter_analysis
  $ProgDir/extract_IG-regions_OR74.py --gff $Gff --fasta $Assembly --prefix $OutDir/OR74A_IG_regions_${Dist} --distance $Dist
  ls $OutDir/OR74A_IG_regions_${Dist}*
done
```


```bash
cd /home/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics
OutDir=analysis/promoters/motifs/N.crassa/OR74A/WC-1/Smith2010
mkdir -p $OutDir

echo "frq vvd cry bli-4 os-4 NCU02264 NCU10063 NCU04021 NCU05594 NCU02800 NCU07541 al-1 NCU03968 NCU11300 al-2" | sed 's/ /\n/g' > $OutDir/gene_list_selected.txt

```


```bash
screen -a
qlogin

ProjDir=$(ls -d /home/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics)
cd $ProjDir

for Promoters in $(ls $ProjDir/analysis/promoters/N.crassa/OR74A/OR74A_IG_regions_*.fa); do
  echo $Promoters
  PromLgth=$(echo $Promoters | rev | cut -f1 -d '_' | rev | sed 's/.fa//g')
  echo $PromLgth  
  OutDir=analysis/promoters/motifs/N.crassa/OR74A/WC-1/Smith2010/selected/${PromLgth}_intergenic_unmasked/meme
  mkdir -p $OutDir

  cat $Promoters | grep '>' | grep -v -f analysis/promoters/motifs/N.crassa/OR74A/WC-1/Smith2010/gene_list_selected.txt | tr -d '>' > $OutDir/gene_list_selected_control.txt
  cat $Promoters | grep '>' | grep 'kal-1' | tr -d '>'>> $OutDir/gene_list_selected_control.txt

  GeneList=$(cat analysis/promoters/motifs/N.crassa/OR74A/WC-1/Smith2010/gene_list_selected.txt)
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
done 2>&1 | tee analysis/promoters/motifs/N.crassa/OR74A/WC-1/Smith2010/selected/intergenic_meme_unmasked.log

cat analysis/promoters/motifs/N.crassa/OR74A/WC-1/Smith2010/selected/intergenic_meme_unmasked.log | grep -e '_IG_regions_' -e 'Motif' | grep -B1 'Motif'| less -S
```

The Smith 2010 motif was only identified in the 5000 bp dataset:

```
/home/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics/analysis/promoters/N.crassa/OR74A/OR74A_IG_regions_3500.fa
        Motif GGTAG MEME-1 regular expression
--
/home/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics/analysis/promoters/N.crassa/OR74A/OR74A_IG_regions_5000.fa
        Motif GATCGAGA MEME-1 regular expression
--
/home/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics/analysis/promoters/N.crassa/OR74A/OR74A_IG_regions_6000.fa
        Motif CBACSACSACS MEME-1 regular expression
```

## Iterations through different promoter region lengths


Firstly, promoter sequences were extracted from the genome using
published gene models

```bash
Assembly=$(ls assembly/external_groups/N.crassa/OR74A/GCF_000182925.2_NC12_genomic.fna)
Gff=$(ls assembly/external_groups/N.crassa/OR74A/GCF_000182925.2_NC12_genomic.gff)

OutDir=analysis/promoters/N.crassa/OR74A
mkdir -p $OutDir
ProgDir=/home/armita/git_repos/emr_repos/scripts/verticillium_pathogenomics/analysis/promoter_analysis
for Dist in $(seq 500 500 6000); do
  echo $Dist
  $ProgDir/extract_promoters_OR74A.py --gff $Gff --fasta $Assembly --prefix $OutDir/OR74A_promoters2_${Dist} --distance $Dist
done

ls $OutDir
```


```bash
cd /home/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics
OutDir=analysis/promoters/motifs/N.crassa/OR74A/WC-1/Smith2010
mkdir -p $OutDir

echo "frq vvd cry bli-4 os-4 NCU02264 NCU10063 NCU04021 NCU05594 NCU02800 NCU07541 al-1 NCU03968 NCU11300 al-2" | sed 's/ /\n/g' > $OutDir/gene_list_selected.txt

```


```bash
screen -a
qlogin

ProjDir=$(ls -d /home/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics)
cd $ProjDir

for Promoters in $(ls $ProjDir/analysis/promoters/N.crassa/OR74A/OR74A_promoters2_*.fa); do
  echo $Promoters
  PromLgth=$(echo $Promoters | rev | cut -f1 -d '_' | rev | sed 's/.fa//g')
  echo $PromLgth
  OutDir=analysis/promoters/motifs/N.crassa/OR74A/WC-1/Smith2010/selected/${PromLgth}_promoters/meme_de
  mkdir -p $OutDir

  GeneList=$(cat analysis/promoters/motifs/N.crassa/OR74A/WC-1/Smith2010/gene_list_selected.txt | sed "s/$/ /g" | tr -d '\n')
  SubPromoters=$OutDir/promoters_OR74A_Smith2010_highpeaks.fa
  for GeneID in $GeneList; do
  cat $Promoters | sed -n "/^>${GeneID}_/,/^>/p" | grep -v "^$" | head -n -1
  done > $SubPromoters

  cat $Promoters | grep '>' | grep -v -f analysis/promoters/motifs/N.crassa/OR74A/WC-1/Smith2010/gene_list_selected.txt | tr -d '>' > $OutDir/gene_list_selected_control.txt
  cat $Promoters | grep '>' | grep 'kal-1' | tr -d '>'>> $OutDir/gene_list_selected_control.txt

  GeneList=$(ls $OutDir/gene_list_selected_control.txt)
  ControlPromoters=$OutDir/promoters_OR74A_Smith2010_-ves.fa
  ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
  $ProgDir/extract_from_fasta.py --fasta $Promoters --headers $GeneList > $ControlPromoters
  ls $ControlPromoters

  fasta-get-markov -m 2 $Promoters $OutDir/background_model.txt

  mkdir $OutDir/meme
  meme -objfun de -neg $ControlPromoters -bfile $OutDir/background_model.txt -dna -mod anr -nmotifs 5 -minw 4 -maxw 12 -revcomp -evt 0.05 -oc $OutDir/meme $SubPromoters

  # meme $SubPromoters -dna -mod anr -nmotifs 5 -minw 6 -maxw 12 -revcomp -evt 0.05 -oc $OutDir/meme

  cat $OutDir/meme/meme.txt | grep -C2 'regular expression'

  #---
  # Running MAST
  #---

  mast $OutDir/meme/meme.xml $SubPromoters -oc $OutDir/mast
  # ls $OutDir/mast_${ListLen}/mast.txt
# done > analysis/promoters/motifs/N.crassa/OR74A/WC-1/Smith2010/meme_unmasked2.log

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
done 2>&1 | tee analysis/promoters/motifs/N.crassa/OR74A/WC-1/Smith2010/selected/promoters_meme_unmasked.log

cat analysis/promoters/motifs/N.crassa/OR74A/WC-1/Smith2010/selected/promoters_meme_unmasked.log | grep -e '_IG_regions_' -e 'Motif' | grep -B1 'Motif'| less -S
```

```
/home/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics/analysis/promoters/N.crassa/OR74A/OR74A_IG_regions_3500.fa
        Motif GGTAG MEME-1 regular expression
--
/home/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics/analysis/promoters/N.crassa/OR74A/OR74A_IG_regions_5000.fa
        Motif GATCGAGA MEME-1 regular expression
--
/home/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics/analysis/promoters/N.crassa/OR74A/OR74A_IG_regions_6000.fa
        Motif CBACSACSACS MEME-1 regular expression
4A/WC-1/Smith2010/selected/promoters_meme_unmasked.log | grep -e '_IG_regions_' -e 'Motif' | grep -B1 'Motif'omoters/motifs/N.crassa/OR7
        Motif CMCSMYCVCC MEME-1 regular expression
        Motif GAGA MEME-2 regular expression
        Motif CHGMTCG MEME-1 regular expression
```
