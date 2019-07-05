Commands used for direct RNAseq analysis Verticillium in response to light conditions.

# Data organisation

dRNAseq data was copied across into the Verticillium project directory
<!--
```bash
ls /oldhpc/home/groups/harrisonlab/project_files/neurospora_crassa/qc_direct_rna
```

Data labelled as 'A' or 'B' was from the same run but the gridion had been restarted. These files were combined:

```bash
  ProjDir=$(ls -d /oldhpc/home/groups/harrisonlab/project_files/neurospora_crassa)
	for RnaSeq in $(ls $ProjDir/qc_direct_rna/*_trim.fastq.gz); do
    Condition=$(basename ${RnaSeq%_trim.fastq.gz})
    Condition=$(echo $Condition | sed "s/.$//g")
    echo $Condition
    OutDir=$ProjDir/qc_rna/minion/dRNAseq/$Condition
    mkdir -p $OutDir
    cat $RnaSeq >> $OutDir/${Condition}_trim_appended.fastq.gz
  done
``` -->


# Pinfish

Transcripts were assembled from dRNAseq data under different light conditions.

This was performed on the new cluster.

```bash
ProjDir=/oldhpc/home/groups/harrisonlab/project_files/neurospora_crassa
cd $ProjDir
```

```bash
  screen -a
  ssh compute03
  ProjDir=$(ls -d /oldhpc/home/groups/harrisonlab/project_files/verticillium_dahliae/clocks)
  WorkDir=/home/armita/prog/pinfish/pipeline-pinfish-analysis
  cd $WorkDir
	for RnaSeq in $(ls $ProjDir/nanopore/dRNAseq/*/*.fastq.gz); do
    Condition=$(basename ${RnaSeq%.fastq.gz})
    echo $Condition
  	RefGenome=$(ls $ProjDir/public_genomes/JR2/Verticillium_dahliaejr2.GCA_000400815.2.dna.toplevel.fa)
	  OutDir=$ProjDir/../pathogenomics/gene_pred/pinfish/V.dahliae/JR2/$Condition
    mkdir -p $OutDir
  	cat config_template.yml | sed "s&genome_fasta:.*&genome_fasta: \"${RefGenome}\"&g" | sed "s&reads_fastq:.*&reads_fastq: \"${RnaSeq}\"&g" | sed "s&workdir_top:.*&workdir_top: \"$OutDir\"&g" > config.yml
    cp config.yml $OutDir
  	snakemake --use-conda -j 20 all
  done
```

 view results
```bash
	ls -lh Workspaces/pipeline-pinfish-analysis/results/
```


# Illumina RNAseq alignment

To assess whether the predicted gene models are correct, RNAseq data was aligned to the JR2 genome.


```bash
for Assembly in $(ls ../clocks/public_genomes/JR2/Verticillium_dahliaejr2.GCA_000400815.2.dna.toplevel.fa); do
Strain=$(echo $Assembly| rev | cut -d '/' -f2 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f3 | rev)
echo "$Organism - $Strain"
for FileF in $(ls qc_rna/paired/*/*/F/*_trim.fq.gz); do
FileR=$(echo $FileF | sed 's&/F/&/R/&g' | sed 's/R1/R2/g')
echo $FileF
echo $FileR
Prefix=$(echo $FileF | rev | cut -f3 -d '/' | rev)
# Timepoint=$(echo $FileF | rev | cut -f2 -d '/' | rev)
Timepoint="timepoint"
#echo "$Timepoint"
OutDir=alignment/star/$Organism/$Strain/$Timepoint/$Prefix
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq
qsub $ProgDir/sub_star.sh $Assembly $FileF $FileR $OutDir
done
done
```

## Run sub_star.sh in a loop

```bash
for FilePath in $(ls -d ../clocks/qc_rna/experiment2/V.dahliae/* | grep 'rep1' | tail -n+2); do
Strain=$(echo $FilePath | rev | cut -f1 -d '/' | rev)
InGenome=$(ls ../clocks/public_genomes/JR2/Verticillium_dahliaejr2.GCA_000400815.2.dna.toplevel.fa)
InGff=$(ls ../clocks/public_genomes/JR2/Verticillium_dahliaejr2.GCA_000400815.2.33.gff3)
InReadF=$(ls ../clocks/qc_rna/experiment2/V.dahliae/$Strain/F/*.fq.gz)
InReadR=$(ls ../clocks/qc_rna/experiment2/V.dahliae/$Strain/R/*.fq.gz)
OutDir=../../../../../../../../data/scratch/armita/verticillium_dahliae/pathogenomics/alignment/STAR/experiment2/V.dahliae/$Strain/
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq
qsub $ProgDir/sub_star.sh $InGenome $InReadF $InReadR $OutDir $InGff
done
```

```bash
for FilePath in $(ls -d ../clocks/qc_rna/experiment1/V.dahliae/* | grep 'rep1'); do
Strain=$(echo $FilePath | rev | cut -f1 -d '/' | rev)
InGenome=$(ls ../clocks/public_genomes/JR2/Verticillium_dahliaejr2.GCA_000400815.2.dna.toplevel.fa)
InGff=$(ls ../clocks/public_genomes/JR2/Verticillium_dahliaejr2.GCA_000400815.2.33.gff3)
InReadF=$(ls $FilePath/F/*.fq.gz)
InReadR=$(ls $FilePath/R/*.fq.gz)
OutDir=../../../../../../../../data/scratch/armita/verticillium_dahliae/pathogenomics/alignment/STAR/experiment1/V.dahliae/$Strain/
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq
qsub $ProgDir/sub_star.sh $InGenome $InReadF $InReadR $OutDir $InGff
done
```

Alignment and the input assembly files were indexed:

```bash
CurDir=$PWD
for File in $(ls alignment/star/public_genomes/JR2/timepoint/*/star_aligment*.bam); do
  cd $(dirname $File)
  bamtools index -in $(basename $File)
  cd $CurDir
done
```

Alignment and the input assembly files were indexed:

```bash
CurDir=$PWD
for File in $(ls /data/scratch/armita/verticillium_dahliae/pathogenomics/alignment/STAR/experiment1/V.dahliae/*/*.bam); do
  cd $(dirname $File)
  bamtools index -in $(basename $File)
  cd $CurDir
done

CurDir=$PWD
for File in $(ls /data/scratch/armita/verticillium_dahliae/pathogenomics/alignment/STAR/experiment2/V.dahliae/*/*.bam); do
  cd $(dirname $File)
  bamtools index -in $(basename $File)
  cd $CurDir
done
```
