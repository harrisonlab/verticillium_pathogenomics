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
RelatedProteins=$(ls $OldProjDir/public_genomes/V.dahliae/JR2/Verticillium_dahliaejr2.GCA_000400815.2.pep.all.fa)
Organism="V.dahliae"
Strain="12008"

OutDir=$OldProjDir/gene_pred/braker/$Organism/${Strain}_UTR
AcceptedHits=$(ls $OldProjDir/alignment/star/V.dahliae/12008/treatment/concatenated/concatenated.bam)
GeneModelName="$Organism"_"$Strain"_braker
CurDir=$PWD

mkdir -p $WorkDir
cd $WorkDir

cp $Assembly assembly.fa
cp $AcceptedHits alignedRNA.bam
cat $RelatedProteins | cut -f1 -d ' ' > related_prots.fa

rm -r /home/armita/prog/augustus/Augustus/config/species/V.dahliae_12008_braker

braker.pl \
  --cores 40 \
  --overwrite \
  --fungus \
  --UTR=on \
  --gff3 \
  --softmasking on \
  --species=$GeneModelName \
  --genome="assembly.fa" \
  # --prot_seq="related_prots.fa" \
  # --prg=gth \
  --bam="alignedRNA.bam"

mkdir -p $OutDir
cp -r braker/* $OutDir/.

rm -r $WorkDir
```


```bash
cd /home/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics
```

Fasta and gff files were extracted from Braker1 output.

```bash
for File in $(ls gene_pred/braker/*/*_UTR/*/augustus.hints.gff3); do
	Strain=$(echo $File | rev | cut -d '/' -f3 | rev | sed 's/_UTR//g')
	Organism=$(echo $File | rev | cut -d '/' -f4 | rev)
	echo "$Organism - $Strain"
	echo "number of genes:"
	cat $File | grep -v '#' | grep -w 'gene' | wc -l
	echo "number of genes with predicted UTRs"
	cat ${File%.gff3}_utr.gff | grep -v '#' | grep -w 'gene' | wc -l

	getAnnoFasta.pl $File
	OutDir=$(dirname $File)
	echo "##gff-version 3" > $OutDir/augustus_extracted.gff
	cat $File | grep -v '#' >> $OutDir/augustus_extracted.gff
done
```

```
  V.dahliae - 12008
  number of genes:
  9602
  number of genes with predicted UTRs
  9196
```


## Antismash secmet prediction with Cassis

Log into the new cluster:

```bash
screen -a
ssh compute03
WorkDir=~/tmp/antismash_Vd
mkdir $WorkDir
cd $WorkDir
OldProjDir=/oldhpc/home/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics
Assembly=$(ls $OldProjDir/repeat_masked/V.*/*/ncbi*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa | grep '12008')
Organism="V.dahliae"
Strain="12008"
Genes=$(ls $OldProjDir/gene_pred/final/$Organism/$Strain/final/final_genes_appended_renamed.gff3)
conda activate emboss
seqret -sequence $Assembly -feature -fformat gff -fopenfile $Genes -osformat genbank -auto -outseq Fv_genes.gbk
conda deactivate

conda activate antismash
antismash \
  --cpus 40 \
  --taxon fungi \
  --transatpks_da \
  --clusterblast \
  --subclusterblast \
  --knownclusterblast \
  --smcogs \
  --inclusive \
  --cassis \
  --borderpredict \
  --full-hmmer \
  --asf \
  --tta \
  --outputfolder Fv_antismash \
  --verbose \
  Fv_genes.gbk
```
