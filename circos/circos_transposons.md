
Script was used to show the positions of effectors, transposons and CAZY gene in 12008 genome.

Sets up variables for Circos

```bash
ProgDir=/home/fanron/git_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos/
Vd12008_genome=repeat_masked/V.dahliae/12008/ncbi_filtered_contigs_repmask/12008_contigs_unmasked.fa
OutDir=analysis/circos/V.dahliae/12008_ticks
mkdir -p $OutDir
```

# Convert the 12008 genome into circos format

$ProgDir/fasta2circos.py --genome $Vd12008_genome --contig_prefix "" > $OutDir/12008_genome.txt

# Make 100kb windows for plots

$ProgDir/fasta2gff_windows.py --genome $Vd12008_genome > $OutDir/12008_100kb_windows.gff

# Convert 12008 MiSeq reads aligning in 100kb windows into coverage stats

```bash
for ReadsBam in $(ls analysis/genome_alignment/bowtie/V.dahliae/12008/vs_12008/12008_contigs_unmasked.fa_aligned_sorted.bam); do
Organism=$(echo $ReadsBam | rev | cut -f4 -d '/' | rev)
Strain=$(echo $ReadsBam | rev | cut -f3 -d '/' | rev)
AlignDir=$(dirname $ReadsBam)
echo "$Organism - $Strain"
bedtools coverage -abam $ReadsBam -b $OutDir/12008_100kb_windows.gff > $AlignDir/"$Strain"_coverage_vs_12008.bed
# Convert coverage bed files into circos format
$ProgDir/coverage_bed2circos.py --bed $AlignDir/"$Strain"_coverage_vs_12008.bed > $OutDir/"$Strain"_coverage_vs_12008_scatterplot.txt 
done
```

# Plot location of 12008 mimps and secreted effectorP genes as a scatterplot
GffTrans=repeat_masked/V.dahliae/12008/ncbi_filtered_contigs_repmask/12008_contigs_transposonmasked.gff
ProgDir=~/git_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
$ProgDir/gff2circos_scatterplot.py --gff $GffTrans --feature similarity --value 1 > $OutDir/12008_Trans1_plot.txt
GffEffP=analysis/effectorP/V.dahliae/12008/V.dahliae_12008_EffectorP_secreted.gff
$ProgDir/gff2circos_scatterplot.py --gff $GffEffP --feature gene --value 0.5 > $OutDir/12008_effectorP_plot.txt
CAZy=gene_pred/CAZY/V.dahliae/12008/12008_CAZY_secreted.gff
$ProgDir/gff2circos_scatterplot.py --gff $CAZy --feature gene --value 0.5 > $OutDir/12008_CAZY_plot.txt
Antismach=analysis/antismash/V.dahliae/12008/12008_secondary_metabolite_regions.gff
$ProgDir/gff2circos_scatterplot.py --gff $Antismach --feature t1pks --value 0.5 > $OutDir/12008_Antismach_plot.txt

circos -conf /home/fanron/git_repos/scripts/verticillium_pathogenomics/circos/circos.conf -outputdir ./$OutDir
```








