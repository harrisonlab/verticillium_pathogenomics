
Commands used in the annotation of the vesca genome for putative resistance
genes as part of resistance gene mapping to V. dahliae


Work was performed in the following directory:

```bash
ProjectDir='/home/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics'
cd $ProjectDir
```

Gene model data was downloaded to the following directory:

```bash
OutDir=gene_pred/downloaded/F.vesca/Hawaii4
mkdir -p $OutDir
CurDir=$PWD
cd $OutDir
wget ftp://ftp.bioinfo.wsu.edu/species/Fragaria_vesca/Fvesca-genome.v1.0/genes/fvesca_v1.0_genemark_hybrid.gff3.gz
wget ftp://ftp.bioinfo.wsu.edu/species/Fragaria_vesca/Fvesca-genome.v1.0/genes/fvesca_v1.0_genemark_hybrid2pseudo.gff3.gz
wget ftp://ftp.bioinfo.wsu.edu/species/Fragaria_vesca/Fvesca-genome.v1.0/genes/fvesca_v1.0_genemark_hybrid.fna.gz
wget ftp://ftp.bioinfo.wsu.edu/species/Fragaria_vesca/Fvesca-genome.v1.0/genes/fvesca_v1.0_genemark_hybrid.annotated.fna.gz
wget ftp://ftp.bioinfo.wsu.edu/species/Fragaria_vesca/Fvesca-genome.v1.0/genes/fvesca_v1.0_genemark_hybrid.faa.gz
wget ftp://ftp.bioinfo.wsu.edu/species/Fragaria_vesca/Fvesca-genome.v1.0/genes/fvesca_v1.0_genemark_hybrid.annotated.faa.gz
wget ftp://ftp.bioinfo.wsu.edu/species/Fragaria_vesca/Fvesca-genome.v1.0/genes/fvesca_v1.0_genemark_abinitio.gff3.gz
wget ftp://ftp.bioinfo.wsu.edu/species/Fragaria_vesca/Fvesca-genome.v1.0/genes/fvesca_v1.0_genemark_abinitio.fna.gz
wget ftp://ftp.bioinfo.wsu.edu/species/Fragaria_vesca/Fvesca-genome.v1.0/genes/fvesca_v1.0_genemark_abinitio.faa.gz
wget ftp://ftp.bioinfo.wsu.edu/species/Fragaria_vesca/Fvesca-genome.v1.1/genes/fvesca_v1.1_genemark_hybrid.gff3.gz
gunzip *.gz
mkdir hybrid
mv *hybrid* hybrid/.
mkdir -p ../Hawaii4_abinitio/abinitio
mv *abinitio* ../Hawaii4_abinitio/abinitio/.
cd $CurDir
```

A combined set of proteins was already available on the cluster:
```bash
ls /home/groups/harrisonlab/project_files/popgen/renseq/strawberry/rgaugury/vesca1.1/lists/fvesca_v1.1_all_proteins.fa
mkdir -p gene_pred/downloaded/F.vesca/Hawaii4_1.1/all_prots
cp /home/groups/harrisonlab/project_files/popgen/renseq/strawberry/rgaugury/vesca1.1/lists/fvesca_v1.1_all_proteins.fa gene_pred/downloaded/F.vesca/Hawaii4_1.1/all_prots/.
cd gene_pred/downloaded/F.vesca/Hawaii4_1.1/all_prots
wget ftp://ftp.bioinfo.wsu.edu/species/Fragaria_vesca/Fvesca-genome.v1.1.a2/genes/Fragaria_vesca_v1.1.a2.gff3.gz
gunzip *.gz
cd $CurDir
```


Helens datafiles were also uploaded to the directory:

```bash
OutDir=analysis/F.vesca/Hawaii4_1.1/helens_files
mkdir -p $OutDir
# Upload files here
```

As such, functional chracterisation was performed on the hybrid gene models of
the 1.0 assembly


```bash
	screen -a
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/interproscan
	# for Genes in $(ls gene_pred/downloaded/F.vesca/Hawaii4/hybrid/fvesca_v1.0_genemark_hybrid.faa); do
	for Genes in $(ls gene_pred/downloaded/F.vesca/Hawaii4_1.1/all_prots/fvesca_v1.1_all_proteins.fa); do
	echo $Genes
	$ProgDir/sub_interproscan.sh $Genes
	done 2>&1 | tee -a interproscan_submisison.log
```

Following interproscan annotation split files were combined using the following
commands:

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/interproscan
for Proteins in $(ls gene_pred/downloaded/F.vesca/Hawaii4_1.1/all_prots/fvesca_v1.1_all_proteins.fa); do
Strain=$(echo $Proteins | rev | cut -d '/' -f3 | rev)
Organism=$(echo $Proteins | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
echo $Strain
InterProRaw=gene_pred/interproscan/$Organism/$Strain/raw
$ProgDir/append_interpro.sh $Proteins $InterProRaw
done
```

proteins containing NB-ARC domains (IPR002182) were identified:

```bash
InterPro=$(ls gene_pred/interproscan/F.vesca/Hawaii4_1.1/Hawaii4_1.1_interproscan.tsv)
OutDir=analysis/F.vesca/Hawaii4_1.1
mkdir -p $OutDir
cat $InterPro | grep 'IPR002182' | cut -f1 | sort | uniq > $OutDir/NB-ARC_genes.txt
cat $OutDir/NB-ARC_genes.txt | sed 's/-1_1//g' | sed 's/hybrid_1/hybrid/g' > $OutDir/NB-ARC_genes_parsed.txt
# GeneGff=$(ls gene_pred/downloaded/F.vesca/Hawaii4_1.1/all_prots/Fragaria_vesca_v1.1.a2.gff3)
# cat $GeneGff | grep -w -f $OutDir/NB-ARC_genes_parsed.txt > $OutDir/NB-ARC_genes.gff
# cat $OutDir/NB-ARC_genes.gff | wc -l
NbsGff=$(ls analysis/F.vesca/Hawaii4_1.1/Fvesca_annotation/vesca1.1.NBS.all.gff3)
cat $NbsGff | grep -w -f $OutDir/NB-ARC_genes_parsed.txt > $OutDir/NB-ARC_genes.gff
cat $OutDir/NB-ARC_genes.gff | wc -l
cat $OutDir/NB-ARC_genes.gff | cut -f9 | cut -f3 -d ';' | sed 's/Parent=//g' > $OutDir/NBS_NB-ARC_genes.txt
cat $InterPro | grep -w -f $OutDir/NBS_NB-ARC_genes.txt > $OutDir/NBS_NB-ARC_interpro.tsv
cat $OutDir/NBS_NB-ARC_interpro.tsv | cut -f1  | sort | uniq | wc -l
```

```
NBS-NB-ARC genes:
290
Extracted headers:
290
Extracted IPR terms for this many extracted headers
196 - somthing has gone wrong
```


Tetratricopeptide repeat-containing domain (IPR013026)
Toll/interleukin-1 receptor homology (TIR) domain (IPR000157)
Virus X resistance protein-like, coiled-coil domain (IPR038005)
WD40-repeat-containing domain (IPR017986)
Powdery mildew resistance protein, RPW8 domain (IPR008808)


## Identifying all domains within 100kb of a SNP

```bash
cat analysis/F.vesca/Hawaii4_1.1/Fvesca_annotation/intersect_output.txt | cut -f19 | cut -f3 -d ';' | grep 'Parent' | sed 's/Parent=//g' | sort | uniq > $OutDir/genes_in_1kb_of_SNP.txt
cat $InterPro | grep -w -f $OutDir/genes_in_1kb_of_SNP.txt > $OutDir/genes_in_1kb_of_SNP_interpro.tsv


cat $OutDir/genes_in_1kb_of_SNP.txt | wc -l
cat gene_pred/downloaded/F.vesca/Hawaii4_1.1/all_prots/fvesca_v1.1_all_proteins.fa | grep '>' | tr -d '>' | grep -w -f $OutDir/genes_in_1kb_of_SNP.txt | wc -l
cat $OutDir/genes_in_1kb_of_SNP_interpro.tsv | cut -f1 | sort | uniq | wc -l


cat $OutDir/genes_in_1kb_of_SNP_interpro.tsv | cut -f1,13 |
sort | uniq
```
