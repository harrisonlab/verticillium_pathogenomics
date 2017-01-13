```bash
  OutDir=analysis/circos/12008_vs_JR2
  mkdir -p $OutDir

  Assembly_12008=repeat_masked/V.dahliae/12008/ncbi_filtered_contigs_repmask/12008_contigs_unmasked.fa
  ProgDir=/home/fanron/git_repos/tools/seq_tools/circos
  $ProgDir/fasta2circos.py --genome $Assembly_12008 --contig_prefix "12008_" > $OutDir/12008_genome.txt

  Assembly_JR2=assembly/merged_canu_spades/V.dahliae/JR2/ensembl/Verticillium_dahliaejr2.GCA_000400815.2.dna.toplevel.fa
  ProgDir=/home/fanron/git_repos/tools/seq_tools/circos
  $ProgDir/fasta2circos.py --genome $Assembly_JR2 --contig_prefix "JR2_" > $OutDir/JR2_genome.txt

  cat $OutDir/12008_genome.txt > $OutDir/12008_JR2_genome.txt
  tac $OutDir/JR2_genome.txt >> $OutDir/12008_JR2_genome.txt

  cat $OutDir/12008_JR2_genome.txt > $OutDir/12008_JR2_genome_edited.txt
```

After an initial round of plotting the order of contigs was manually changed using nano

```bash
nano $OutDir/12008_JR2_genome_edited.txt
```

links were made between single copy orthologous genes. These links were ordered by
JR2 chromosome order to aid the re-ordering of contigs to plot in the
$OutDir/12008_JR2_genome_edited.txt file.

```bash
  ProgDir=/home/fanron/git_repos/scripts/verticillium_pathogenomics/circos/12008_vs_reference
  $ProgDir/orthology2circos_ribbons.py \
  --orthology analysis/orthology/orthomcl/All_Strains/All_Strains_orthogroups.txt \
  --name1 12008 \
  --gff1 gene_pred/final_genes/V.dahliae/12008/final/final_genes_appended.gff3 \
  --name2 JR2 \
  --gff2 assembly/merged_canu_spades/V.dahliae/JR2/ensembl/Verticillium_dahliaejr2.GCA_000400815.2.33.gff3 \
  | sort -k4,5 -V > $OutDir/12008_JR2_links.txt
  cat $OutDir/12008_JR2_links.txt > $OutDir/12008_JR2_links_edited.txt
  ProgDir=/home/fanron/git_repos/scripts/verticillium_pathogenomics/circos/12008_vs_reference
  circos -conf $ProgDir/12008_vs_JR2_circos.conf -outputdir $OutDir
  mv $OutDir/circos.png $OutDir/12008_JR2_circos.png
  mv $OutDir/circos.svg $OutDir/12008_JR2_circos.svg
```

a file showing contig orders was made using:

```bash
  cat $OutDir/12008_JR2_links_edited.txt | cut -f1 | uniq > $OutDir/12008_contig_order.txt
```


```bash
OutDir=analysis/circos/F.oxysporum_fsp_cepae/Fus2_FoL
cat $OutDir/12008_JR2_genome_edited2.txt | grep -v '4287' > $OutDir/12008_JR2_genome_final.txt
mkdir -p $OutDir/by_FoC_chr
for Num in $(seq 1 22); do
  Chr="contig_"$Num"_pilon"
  echo "$Chr"
  OrthologyTxt=analysis/orthology/orthomcl/FoC_vs_Fo_vs_FoL_publication_ncbi/FoC_vs_Fo_vs_FoL_publication_ncbi_orthogroups.txt
  ProgDir=/home/armita/git_repos/emr_repos/scripts/verticillium_pathogenomics/circos/12008_vs_reference
  $ProgDir/orthology2ribbons_internal.py \
  --chr1 $Chr \
  --orthology $OrthologyTxt \
  --name1 Fus2 \
  --gff1 gene_pred/final_genes/F.oxysporum_fsp_cepae/Fus2_canu_new/final/final_genes_appended.gff3 \
  | sort | uniq \
  > $OutDir/12008_JR2_links_edited.txt
  ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  circos -conf $ProgDir/Fus2/Fus2_FoL/12008_JR2_circos.conf -outputdir $OutDir
  mv $OutDir/circos.png $OutDir/by_FoC_chr/12008_JR2_LS_"$Chr"_circos.png
  mv $OutDir/circos.svg $OutDir/by_FoC_chr/12008_JR2_LS_"$Chr"_circos.svg
done
```

The frequency of gene duplications within and between chromosomes was investigated:

```bash
OutDir=analysis/circos/F.oxysporum_fsp_cepae/Fus2_FoL
for Num in $(seq 1 22); do
Chr="contig_"$Num"_pilon"
ChrList="$ChrList $Chr"
done
echo "$ChrList"
OrthologyTxt=analysis/orthology/orthomcl/FoC_vs_Fo_vs_FoL_publication_ncbi/FoC_vs_Fo_vs_FoL_publication_ncbi_orthogroups.txt
ProgDir=/home/armita/git_repos/emr_repos/scripts/verticillium_pathogenomics/circos/12008_vs_reference
$ProgDir/orthology2ribbons_internal.py \
--chr1 $ChrList \
--orthology $OrthologyTxt \
--name1 Fus2 \
--gff1 gene_pred/final_genes/F.oxysporum_fsp_cepae/Fus2_canu_new/final/final_genes_appended.gff3 \
| sort | uniq > $OutDir/Fus2_all_Fus2_links.txt
cat $OutDir/Fus2_all_Fus2_links.txt | cut -f1,4 | sort | uniq -c > $OutDir/Fus2_all_Fus2_links_occurence.txt

```
