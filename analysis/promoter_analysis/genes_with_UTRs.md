# Genes with UTRs

Gene models were repredicted using an updated version of Braker
on the new cluster that allows UTR prediction.

This location is in the /data2 area of the old cluster:

```bash
  ProjDir=/oldhpc/data/scratch/armita/verticillium_dahliae/pathogenomics
  mkdir -p $ProjDir
  cd $ProjDir
```

Gene prediction was run from a screen session with a ssh connection to a worker node:

```bash
screen -a
ssh compute02
cd /oldhpc/data/scratch/armita/verticillium_dahliae/pathogenomics


WorkDir=$HOME/tmp/braker_Vd
OldProjDir=/oldhpc/home/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics

Assembly=$(ls $OldProjDir/repeat_masked/V.*/*/ncbi*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa | grep '12008')
Organism="V.dahliae"
Strain="12008"

OutDir=/oldhpc/data/scratch/armita/fusarium_venenatum/fusarium_venenatum/gene_pred/braker/F.venenatum/WT_braker_UTR
AcceptedHits=$(ls $OldProjDir/alignment/star/V.dahliae/12008/treatment/concatenated/concatenated.bam)
GeneModelName="$Organism"_"$Strain"_braker
CurDir=$PWD

mkdir -p $WorkDir
cd $WorkDir

cp $Assembly assembly.fa
cp $AcceptedHits alignedRNA.bam

braker.pl \
  --cores 40 \
  --overwrite \
  --fungus \
  --UTR=on \
  --gff3 \
  --softmasking on \
  --species=$GeneModelName \
  --genome="assembly.fa" \
  --bam="alignedRNA.bam"

mkdir -p $CurDir/$OutDir
cp -r braker/* $CurDir/$OutDir/.

rm -r $WorkDir

```
