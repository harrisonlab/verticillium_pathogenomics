
For all 5 isolates, after bowtie alignment:

Sets up variables for Circos

```bash
ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos/
Vd12008_genome=repeat_masked/V.dahliae/12008/ncbi_filtered_contigs_repmask/12008_contigs_unmasked.fa
OutDir=analysis/circos/V.dahliae/12008
mkdir -p $OutDir
```

# Convert the Fus2 genome into circos format

$ProgDir/fasta2circos.py --genome $Vd12008_genome --contig_prefix "" > $OutDir/12008_merged1_assembly_genome.txt

# Make 100kb windows for plots

$ProgDir/fasta2gff_windows.py --genome $Vd12008_genome > $OutDir/12008_merged1_assembly_100kb_windows.gff

# Convert FoC MiSeq reads aligning in 100kb windows into coverage stats

```bash
for ReadsBam in $(ls analysis/genome_alignment/bowtie/V.dahliae/51/vs_12008_unmasked/12008_contigs_unmasked.fa_aligned_sorted.bam); do
Organism=$(echo $ReadsBam | rev | cut -f4 -d '/' | rev)
Strain=$(echo $ReadsBam | rev | cut -f3 -d '/' | rev)
AlignDir=$(dirname $ReadsBam)
echo "$Organism - $Strain"
bedtools coverage -abam $ReadsBam -b $OutDir/12008_merged1_assembly_100kb_windows.gff > $AlignDir/"$Strain"_coverage_vs1_12008.bed
# Convert coverage bed files into circos format
$ProgDir/coverage_bed2circos.py --bed $AlignDir/"$Strain"_coverage_vs1_12008.bed > $OutDir/"$Strain"_coverage_vs1_12008_scatterplot.txt 
done
```

# Plot location of Fus2 mimps and secreted effectorP genes as a scatterplot
GffTrans=repeat_masked/V.dahliae/12008/ncbi_filtered_contigs_repmask/12008_contigs_transposonmasked.gff
ProgDir=~/git_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
$ProgDir/gff2circos_scatterplot.py --gff $GffTrans --feature similarity --value 1 > $OutDir/12008_Trans1_plot.txt
GffEffP=analysis/effectorP/V.dahliae/12008/*_EffectorP_secreted.gff
$ProgDir/gff2circos_scatterplot.py --gff $GffEffP --feature gene --value 0.5 > $OutDir/R0905_effectorP_plot.txt
circos -conf /home/fanron/git_repos/scripts/verticillium_pathogenomics/circos/circos.conf -outputdir ./$OutDir
```








