(correct) Circos plots were generated using unmasked genome:

```bash
ProgDir=/home/fanron/git_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
Vd12008_genome=repeat_masked/*/*/*/*_contigs_unmasked.fa
OutDir=analysis/circos/V.dahliae/12008_merged1_assembly
mkdir -p "$OutDir"


# Convert the Fus2 genome into circos format
$ProgDir/fasta2circos.py --genome $Vd12008_genome --contig_prefix "" > $OutDir/12008_merged1_assembly_genome.txt

# Make 100kb windows for plots
$ProgDir/fasta2gff_windows.py --genome $Vd12008_genome --size 20 > $OutDir/12008_merged1_assembly_20kb_windows.gff

# Convert FoC MiSeq reads aligning in 100kb windows into coverage stats
```bash
for ReadsBam in $(ls analysis/genome_alignment/bowtie/*/*/*/12008_contigs_unmasked.fa_aligned_sorted.bam); do
Organism=$(echo $ReadsBam | rev | cut -f4 -d '/' | rev)
Strain=$(echo $ReadsBam | rev | cut -f3 -d '/' | rev)
AlignDir=$(dirname $ReadsBam)
echo "$Organism - $Strain"
bedtools coverage -abam $ReadsBam -b $OutDir/12008_merged1_assembly_20kb_windows.gff > $AlignDir/"$Strain"_coverage_vs1_12008.bed
# Convert coverage bed files into circos format
$ProgDir/coverage_bed2circos.py --bed $AlignDir/"$Strain"_coverage_vs1_12008.bed > $OutDir/"$Strain"_coverage_vs1_12008_scatterplot.txt 
done
```

# Plot location of Fus2 mimps and secreted effectorP genes as a scatterplot
	GffTrans=repeat_masked/V.dahliae/12008/filtered_contigs_repmask/12008_contigs_transposonmasked.gff
	ProgDir=~/git_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
	$ProgDir/gff2circos_scatterplot.py --gff $GffTrans --feature similarity --value 1 > $OutDir/12008_Trans1_plot.txt
	#$ProgDir/gff2circos_scatterplot.py --gff $GffTrans --feature similarity --value 0.5 > $OutDir/12008_Trans2_plot.txt
	#GffEffP=analysis/effectorP/N.ditissima/R0905_merged_assembly/N.ditissima_R0905_merged_assembly_EffectorP_secreted.gff
	#ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
	#$ProgDir/gff2circos_scatterplot.py --gff $GffEffP --feature gene --value 0.5 > $OutDir/R0905_effectorP_plot.txt
	circos -conf /home/fanron/git_repos/scripts/verticillium_pathogenomics/circos/circos.conf -outputdir ./$OutDir
```


(Incorrect) Circos plots were generated using softmasked genome:

```bash
ProgDir=/home/fanron/git_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
Vd12008_genome=repeat_masked/*/*/*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa
OutDir=analysis/circos/V.dahliae/12008_merged_assembly
mkdir -p "$OutDir"

# Convert the Fus2 genome into circos format
$ProgDir/fasta2circos.py --genome $Vd12008_genome --contig_prefix "" > $OutDir/12008_merged_assembly_genome.txt

# Make 100kb windows for plots
$ProgDir/fasta2gff_windows.py --genome $Vd12008_genome --size 20 > $OutDir/12008_merged_assembly_20kb_windows.gff

# Convert FoC MiSeq reads aligning in 100kb windows into coverage stats
for ReadsBam in $(ls analysis/genome_alignment/bowtie/*/*/vs_12008softmaskd/12008_contigs_softmasked_repeatmasker_TPSI_appended.fa_aligned.bam); do
Organism=$(echo $ReadsBam | rev | cut -f4 -d '/' | rev)
Strain=$(echo $ReadsBam | rev | cut -f3 -d '/' | rev)
AlignDir=$(dirname $ReadsBam)
echo "$Organism - $Strain"
bedtools coverage -abam $ReadsBam -b $OutDir/12008_merged_assembly_20kb_windows.gff > $AlignDir/"$Strain"_coverage_vs_12008.bed
# Convert coverage bed files into circos format
$ProgDir/coverage_bed2circos.py --bed $AlignDir/"$Strain"_coverage_vs_12008.bed > $OutDir/"$Strain"_coverage_vs_12008_scatterplot.txt
done

# Plot location of Fus2 mimps and secreted effectorP genes as a scatterplot
GffTrans=repeat_masked/V.dahliae/12008/filtered_contigs_repmask/12008_contigs_transposonmasked.gff
ProgDir=~/git_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
$ProgDir/gff2circos_scatterplot.py --gff $GffTrans --feature gene --value 1 > $OutDir/12008_Trans_plot.txt
GffEffP=analysis/effectorP/N.ditissima/R0905_merged_assembly/N.ditissima_R0905_merged_assembly_EffectorP_secreted.gff
ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
$ProgDir/gff2circos_scatterplot.py --gff $GffEffP --feature gene --value 0.5 > $OutDir/R0905_effectorP_plot.txt

*****
circos -conf /home/fanron/git_repos/scripts/verticillium_pathogenomics/circos/circos.conf -outputdir ./$OutDir
```










