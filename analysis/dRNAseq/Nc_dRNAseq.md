Commands used for direct RNAseq analysis Neuropspora in response to light conditions.

# Data organisation


dRNAseq data has copied across into the Neurospora project directory and trimmed by Suresh

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
```


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
  ProjDir=$(ls -d /oldhpc/home/groups/harrisonlab/project_files/neurospora_crassa)
  WorkDir=/home/armita/prog/pinfish/pipeline-pinfish-analysis
  cd $WorkDir
	for RnaSeq in $(ls $ProjDir/qc_rna/minion/dRNAseq/*/*_trim_appended.fastq.gz); do
	# for RnaSeq in $(ls $ProjDir/../verticillium_dahliae/clocks/nanopore/dRNAseq/*/*.fastq.gz); do
    # Condition=$(echo $RnaSeq | rev | cut -f2 -d '/' | rev)
    # Condition=$(basename ${RnaSeq%.fastq.gz})
    Condition=$(basename ${RnaSeq%_trim_appended.fastq.gz})
    echo $Condition
  	RefGenome=$(ls $ProjDir/Neurospora_crassa.NC12.dna_rm.toplevel.fa)
	  OutDir=$ProjDir/gene_pred/pinfish/N.crassa/OR74A/$Condition
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
