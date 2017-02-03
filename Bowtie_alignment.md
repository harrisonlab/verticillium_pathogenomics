##Alignment of raw reads vs the 12008 genome

```bash
Reference=$(ls repeat_masked/V.dahliae/12008/ncbi*/*_contigs_unmasked.fa)
for StrainPath in $(ls -d qc_dna/V.dahliae/paired/*/*)
do
echo $StrainPath
Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
echo "$Organism - $Strain"
F_Read=$(ls $StrainPath/F/*_trim.fq.gz | head -n1 | tail -n1)
R_Read=$(ls $StrainPath/R/*_trim.fq.gz | head -n1 | tail -n1)
echo $F_Read
echo $R_Read
OutDir=analysis/genome_alignment/bowtie/$Organism/$Strain/vs_12008unmasked
ProgDir=/home/fanron/git_repos/tools/seq_tools/genome_alignment
qsub $ProgDir/bowtie/sub_bowtie.sh $Reference $F_Read $R_Read $OutDir $Strain
done
```




##Alignment of raw reads vs the JR2 genome
Reference=$(ls assembly/merged_canu_spades/V.dahliae/JR2/ensembl/Verticillium_dahliaejr2.GCA_000400815.2.dna.toplevel.fa)
for StrainPath in $(ls -d qc_dna/paired/V.dahliae/*/* | grep -e '12008' -e '51')
do
echo $StrainPath
Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
echo "$Organism - $Strain"
F_Read=$(ls $StrainPath/F/*_trim.fq.gz | head -n1 | tail -n1)
R_Read=$(ls $StrainPath/R/*_trim.fq.gz | head -n1 | tail -n1)
echo $F_Read
echo $R_Read
OutDir=analysis/genome_alignment/bowtie/$Organism/$Strain/vs_JR2
ProgDir=/home/fanron/git_repos/tools/seq_tools/genome_alignment
qsub $ProgDir/bowtie/sub_bowtie.sh $Reference $F_Read $R_Read $OutDir $Strain
done
```

##Alignment of 12008 Miseq reads vs the 12008 genome to see the assembly status.

Reference=$(ls repeat_masked/V.dahliae/12008/ncbi*/*_contigs_unmasked.fa)
for StrainPath in $(ls -d qc_dna/paired/V.dahliae/12008)
do
echo $StrainPath
Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
echo "$Organism - $Strain"
F_Read=$(ls $StrainPath/F/*_trim.fq.gz | head -n1 | tail -n1)
R_Read=$(ls $StrainPath/R/*_trim.fq.gz | head -n1 | tail -n1)
echo $F_Read
echo $R_Read
OutDir=analysis/genome_alignment/bowtie/$Organism/$Strain/vs_12008unmasked
ProgDir=/home/fanron/git_repos/tools/seq_tools/genome_alignment
qsub $ProgDir/bowtie/sub_bowtie.sh $Reference $F_Read $R_Read $OutDir $Strain
done