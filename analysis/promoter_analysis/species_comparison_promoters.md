# Specieis comparison motif identification

Commands used to identify regulatory motifs from comparison of promoter regions
of given genes between species.

Where motifs are identified in regualtory element of a given species then this
can be searched in wider gene sets that have a similar expression pattern.

# Target Species

High confidence:
```
Neurospora crassa - clock
Botrytis cinerea - clock
Magnaporthe oryzae - clock
Pyronema confluens - clock
```

Data was downloaded to the cluster:

```bash
ProjDir=/home/groups/harrisonlab/project_files/verticillium_dahliae/pathogenomics
# N. crassa OR74A
Organism="N.crassa"
Strain="OR74A"
OutDir=assembly/external_groups/$Organism/$Strain
mkdir -p $OutDir
Assembly=$(ls assembly/external_groups/N.crassa/OR74A/GCF_000182925.2_NC12_genomic.fna)
Gff=$(ls assembly/external_groups/N.crassa/OR74A/GCF_000182925.2_NC12_genomic.gff)

# B. cinerea
Organism="B.cinerea"
Strain="BO5_10"
OutDir=assembly/external_groups/$Organism/$Strain
mkdir -p $OutDir
cd $OutDir
wget ftp://ftp.ensemblgenomes.org/pub/fungi/release-43/fasta/botrytis_cinerea/dna/*.fa.gz
wget ftp://ftp.ensemblgenomes.org/pub/fungi/release-43/gff3/botrytis_cinerea/*.gff3.gz
gunzip *.gz
cd $ProjDir
Assembly=$(ls assembly/external_groups/B.cinerea/BO5_10/Botrytis_cinerea.ASM83294v1.dna.toplevel.fa)
Gff=$(ls assembly/external_groups/B.cinerea/BO5_10/Botrytis_cinerea.ASM83294v1.43.gff3)

# M. oryzae
Organism="M.oryzae"
Strain="MG08"
OutDir=assembly/external_groups/$Organism/$Strain
mkdir -p $OutDir
cd $OutDir
wget ftp://ftp.ensemblgenomes.org/pub/fungi/release-43/fasta/magnaporthe_oryzae/dna/*.fa.gz
wget ftp://ftp.ensemblgenomes.org/pub/fungi/release-43/gff3/magnaporthe_oryzae/*.gff3.gz
gunzip *.gz
cd $ProjDir
Assembly=$(ls assembly/external_groups/M.oryzae/MG08/Magnaporthe_oryzae.MG8.dna.toplevel.fa)
Gff=$(ls assembly/external_groups/M.oryzae/MG08/Magnaporthe_oryzae.MG8.43.gff3)

# Pyronema confluens
Organism="P.omphalodes"
Strain="CBS100304"
OutDir=assembly/external_groups/$Organism/$Strain
mkdir -p $OutDir
# Downloaded from https://genome.jgi.doe.gov/portal/pages/dynamicOrganismDownload.jsf?organism=Pyrom1

gunzip $OutDir/*.gz
cd $ProjDir
Assembly=$(ls assembly/external_groups/P.omphalodes/CBS100304/Pyrom1_AssemblyScaffolds.fasta)
Gff=$(ls assembly/external_groups/P.omphalodes/CBS100304/Pyrom1_GeneCatalog_genes_20181022.gff)
```
Low confidence:

```
Aureobasidium pullulans - conidial banding
Zymoseptoria trici - conidial banding
Sclerotinia fructigena - conidial banding - no gene models
Sordaria fimicola - conidial banding - no genome
Trichodera pleuroticola - conidial banding - no genome
```

Low confidence:

```bash
# Aureobasidium pullulans - conidial banding
Organism="A.pullulans"
Strain="BO5_10"
OutDir=assembly/external_groups/$Organism/$Strain
mkdir -p $OutDir
cd $OutDir
wget ftp://ftp.ensemblgenomes.org/pub/fungi/release-43/fasta/fungi_ascomycota1_collection/aureobasidium_pullulans_exf_150_gca_000721785/dna/*.fa.gz
wget ftp://ftp.ensemblgenomes.org/pub/fungi/release-43/gff3/fungi_ascomycota1_collection/aureobasidium_pullulans_exf_150_gca_000721785/*.gff3.gz
gunzip *.gz
cd $ProjDir
Assembly=$(ls assembly/external_groups/A.pullulans/BO5_10/Aureobasidium_pullulans_exf_150_gca_000721785.Aureobasidium_pullulans_var._pullulans_EXF-150_assembly_version_1.0.dna.toplevel.fa)
Gff=$(ls assembly/external_groups/A.pullulans/BO5_10/Aureobasidium_pullulans_exf_150_gca_000721785.Aureobasidium_pullulans_var._pullulans_EXF-150_assembly_version_1.0.43.gff3)

# Zymoseptoria trici - conidial banding
Organism="Z.tritici"
Strain="MYCGR_v2"
OutDir=assembly/external_groups/$Organism/$Strain
mkdir -p $OutDir
cd $OutDir
wget ftp://ftp.ensemblgenomes.org/pub/fungi/release-44/fasta/zymoseptoria_tritici/dna/*.fa.gz
wget ftp://ftp.ensemblgenomes.org/pub/fungi/release-44/gff3/zymoseptoria_tritici/*.gff3.gz
gunzip *.gz
cd $ProjDir
Assembly=$(ls assembly/external_groups/Z.tritici/MYCGR_v2/Zymoseptoria_tritici.MG2.dna.toplevel.fa)
Gff=$(ls assembly/external_groups/Z.tritici/MYCGR_v2/Zymoseptoria_tritici.MG2.44.gff3)

# # Sclerotinia (Monilina) fructigena - conidial banding
# Organism="M.fructigena"
# Strain="Mfrg269"
# OutDir=assembly/external_groups/$Organism/$Strain
# mkdir -p $OutDir
# cd $OutDir
# wget https://sra-download.ncbi.nlm.nih.gov/traces/wgs03/wgs_aux/QK/RW/QKRW01/QKRW01.1.fsa_nt.gz
# # Gene models have been predicted but are on ncbi as proteins rather than annotations of the genome

# Sordaria fimicola - conidial banding
# no genome - Sordaria macrosporia genome avalailable
# https://www.ncbi.nlm.nih.gov/genome/?term=Sordaria

# Trichoderma pleuroticola - conidial banding
# no genome available
# other trichoderma genomes available


```

## Extract promoters


### Clock organisms

Firstly, promoter sequences were extracted from the genome using
published gene models

```bash
# Assembly=$(ls ../clocks/public_genomes/JR2/Verticillium_dahliaejr2.GCA_000400815.2.dna.toplevel_parsed.fa)
Assembly=$(ls assembly/external_groups/N.crassa/OR74A/GCF_000182925.2_NC12_genomic.fna)
cat $Assembly | cut -f1 -d ' '> ${Assembly%.fna}_parsed.fa

Gff=$(ls assembly/external_groups/N.crassa/OR74A/GCF_000182925.2_NC12_genomic.gff)
Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
OutDir=analysis/promoters/$Organism/$Strain
mkdir -p $OutDir
Dist=3000
# ProgDir=/home/armita/git_repos/emr_repos/scripts/verticillium_pathogenomics/analysis/promoter_analysis
# $ProgDir/extract_promoters_ensembl.py --fasta ${Assembly%.fa}_parsed.fa --gff $Gff --prefix $OutDir/${Organism}_promoters_${Dist} --distance $Dist

ls $OutDir
```

```bash
# Assembly=$(ls ../clocks/public_genomes/JR2/Verticillium_dahliaejr2.GCA_000400815.2.dna.toplevel_parsed.fa)
Assembly=$(ls assembly/external_groups/B.cinerea/BO5_10/Botrytis_cinerea.ASM83294v1.dna.toplevel.fa)
cat $Assembly | cut -f1 > ${Assembly%.fa}_parsed.fa

Gff=$(ls assembly/external_groups/B.cinerea/BO5_10/Botrytis_cinerea.ASM83294v1.43.gff3)
Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
OutDir=analysis/promoters/$Organism/$Strain
mkdir -p $OutDir
Dist=3000
ProgDir=/home/armita/git_repos/emr_repos/scripts/verticillium_pathogenomics/analysis/promoter_analysis
$ProgDir/extract_promoters_ensembl.py --fasta ${Assembly%.fa}_parsed.fa --gff $Gff --prefix $OutDir/${Organism}_promoters_${Dist} --distance $Dist

ls $OutDir
```

```bash
# Assembly=$(ls ../clocks/public_genomes/JR2/Verticillium_dahliaejr2.GCA_000400815.2.dna.toplevel_parsed.fa)
Assembly=$(ls assembly/external_groups/B.cinerea/BO5_10/Botrytis_cinerea.ASM83294v1.dna.toplevel.fa)
cat $Assembly | cut -f1 > ${Assembly%.fa}_parsed.fa

Gff=$(ls assembly/external_groups/B.cinerea/BO5_10/Botrytis_cinerea.ASM83294v1.43.gff3)
Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
OutDir=analysis/promoters/$Organism/$Strain
mkdir -p $OutDir
Dist=3000
ProgDir=/home/armita/git_repos/emr_repos/scripts/verticillium_pathogenomics/analysis/promoter_analysis
$ProgDir/extract_promoters_ensembl.py --fasta ${Assembly%.fa}_parsed.fa --gff $Gff --prefix $OutDir/${Organism}_promoters_${Dist} --distance $Dist

ls $OutDir
```


## Frq
