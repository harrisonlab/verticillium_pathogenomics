
The final number of genes per isolate was observed using:
```bash
for DirPath in $(ls -d assembly/merged_canu_spades/V.dahliae/*/ensembl); do
echo $DirPath
cat $DirPath/*.pep.all.fa | grep '>' | wc -l
echo ""
done
```
assembly/merged_canu_spades/V.dahliae/JR2/ensembl
11424

assembly/merged_canu_spades/V.dahliae/Ls17/ensembl
10535

## Effector genes

Putative pathogenicity and effector related genes were identified within Braker
gene models using a number of approaches:

 * A) From Augustus gene models - Identifying secreted proteins
 * B) From Augustus gene models - Effector identification using EffectorP


## A) From Augustus gene models - Identifying secreted proteins

 Required programs:
  * SignalP-4.1
  * TMHMM

 Proteins that were predicted to contain signal peptides were identified using
 the following commands:

 ```bash
  SplitfileDir=/home/fanron/git_repos/tools/seq_tools/feature_annotation/signal_peptides
  ProgDir=/home/fanron/git_repos/tools/seq_tools/feature_annotation/signal_peptides
  CurPath=$PWD
  for Proteome in $(ls assembly/merged_canu_spades/V.dahliae/JR2/ensembl/*.pep.all.fa); do
    Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
    SplitDir=gene_pred/final_genes_split/$Organism/$Strain
    mkdir -p $SplitDir
    BaseName="$Organism""_$Strain"_final_preds
    $SplitfileDir/splitfile_500.py --inp_fasta $Proteome --out_dir $SplitDir --out_base $BaseName
    for File in $(ls $SplitDir/*_final_preds_*); do
      Jobs=$(qstat | grep 'pred_sigP' | grep 'qw' | wc -l)
      while [ $Jobs -gt 1 ]; do
        sleep 10
        printf "."
        Jobs=$(qstat | grep 'pred_sigP' | grep 'qw' | wc -l)
      done
      printf "\n"
      echo $File
      qsub $ProgDir/pred_sigP.sh $File signalp-4.1
    done
  done
 ```

 ```bash
  SplitfileDir=/home/fanron/git_repos/tools/seq_tools/feature_annotation/signal_peptides
  ProgDir=/home/fanron/git_repos/tools/seq_tools/feature_annotation/signal_peptides
  CurPath=$PWD
  for Proteome in $(ls assembly/merged_canu_spades/V.alfa*/*/ensembl/*.pep.all.fa); do
    Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
    SplitDir=gene_pred/final_genes_split/$Organism/$Strain
    mkdir -p $SplitDir
    BaseName="$Organism""_$Strain"_final_preds
    $SplitfileDir/splitfile_500.py --inp_fasta $Proteome --out_dir $SplitDir --out_base $BaseName
    for File in $(ls $SplitDir/*_final_preds_*); do
      Jobs=$(qstat | grep 'pred_sigP' | grep 'qw' | wc -l)
      while [ $Jobs -gt 1 ]; do
        sleep 10
        printf "."
        Jobs=$(qstat | grep 'pred_sigP' | grep 'qw' | wc -l)
      done
      printf "\n"
      echo $File
      qsub $ProgDir/pred_sigP.sh $File signalp-4.1
    done
  done
 ```
 The batch files of predicted secreted proteins needed to be combined into a
 single file for each strain. This was done with the following commands:
 
 ```bash
  for SplitDir in $(ls -d gene_pred/final_genes_split/*/*); do
    Strain=$(echo $SplitDir | rev |cut -d '/' -f1 | rev)
    Organism=$(echo $SplitDir | rev |cut -d '/' -f2 | rev)
    InStringAA=''
    InStringNeg=''
    InStringTab=''
    InStringTxt=''
    SigpDir=final_genes_signalp-4.1
    for GRP in $(ls -l $SplitDir/*_final_preds_*.fa | rev | cut -d '_' -f1 | rev | sort -n); do
      InStringAA="$InStringAA gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_final_preds_$GRP""_sp.aa";
      InStringNeg="$InStringNeg gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_final_preds_$GRP""_sp_neg.aa";
      InStringTab="$InStringTab gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_final_preds_$GRP""_sp.tab";
      InStringTxt="$InStringTxt gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_final_preds_$GRP""_sp.txt";
    done
    cat $InStringAA > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_final_sp.aa
    cat $InStringNeg > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_final_neg_sp.aa
    tail -n +2 -q $InStringTab > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_final_sp.tab
    cat $InStringTxt > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_final_sp.txt
  done
 ```
 Some proteins that are incorporated into the cell membrane require secretion.
 Therefore proteins with a transmembrane domain are not likely to represent
 cytoplasmic or apoplastic effectors.

 Proteins containing a transmembrane domain were identified:

 ```bash
  for Proteome in $(ls assembly/merged_canu_spades/*/*/ensembl/*.pep.all.fa); do
    Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
    ProgDir=/home/fanron/git_repos/tools/seq_tools/feature_annotation/transmembrane_helices
    qsub $ProgDir/submit_TMHMM.sh $Proteome
  done
 ```

 Those proteins with transmembrane domains were removed from lists of Signal peptide containing proteins

```bash
  for File in $(ls gene_pred/trans_mem/*/*/*_TM_genes_neg.txt); do
    Strain=$(echo $File | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $File | rev | cut -f3 -d '/' | rev)
    echo "$Organism - $Strain"
    TmHeaders=$(echo "$File" | sed 's/neg.txt/neg_headers.txt/g')
    cat $File | cut -f1 > $TmHeaders
    SigP=$(ls gene_pred/final_genes_signalp-4.1/$Organism/$Strain/*_final_sp.aa)
    OutDir=$(dirname $SigP)
    ProgDir=/home/fanron/git_repos/tools/gene_prediction/ORF_finder
    $ProgDir/extract_from_fasta.py --fasta $SigP --headers $TmHeaders > $OutDir/"$Strain"_final_sp_no_trans_mem.aa
    cat $OutDir/"$Strain"_final_sp_no_trans_mem.aa | grep '>' | wc -l
  done
```

V.alfafae - VaMs102
866
V.dahliae - 12008
943
V.dahliae - JR2
867
V.dahliae - Ls17
908

## B) From Augustus gene models - Effector identification using EffectorP

Required programs:
 * EffectorP.py

```bash
  for Proteome in $(ls assembly/merged_canu_spades/*/JR2/ensembl/*.pep.all.fa); do
    Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
    BaseName="$Organism"_"$Strain"_EffectorP
    OutDir=analysis/effectorP/$Organism/$Strain
    ProgDir=~/git_repos/tools/seq_tools/feature_annotation/fungal_effectors
    qsub $ProgDir/pred_effectorP.sh $Proteome $BaseName $OutDir
  done
```

********
Those genes that were predicted as secreted and tested positive by effectorP were identified:

```bash
    for File in $(ls analysis/effectorP/V.*/VaMs102/*_EffectorP.txt); do
    Strain=$(echo $File | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $File | rev | cut -f3 -d '/' | rev)
    echo "$Organism - $Strain"
    Headers=$(echo "$File" | sed 's/_EffectorP.txt/_EffectorP_headers.txt/g')
    cat $File | grep 'Effector' | cut -f1 -d ' ' > $Headers
    Secretome=$(ls gene_pred/final_genes_signalp-4.1/V.*/VaMs102/*_final_sp_no_trans_mem.aa)
    OutFile=$(echo "$File" | sed 's/_EffectorP.txt/_EffectorP_secreted.aa/g')
    ProgDir=/home/fanron/git_repos/tools/gene_prediction/ORF_finder
    $ProgDir/extract_from_fasta.py --fasta $Secretome --headers $Headers > $OutFile
    OutFileHeaders=$(echo "$File" | sed 's/_EffectorP.txt/_EffectorP_secreted_headers.txt/g')
    cat $OutFile | grep '>' | tr -d '>' > $OutFileHeaders
    cat $OutFileHeaders | wc -l
    done
      #Gff=$(ls gene_pred/codingquary1/V.dahliae/12008/*/final_genes_appended.gff3)
    #EffectorP_Gff=$(echo "$File" | sed 's/_EffectorP.txt/_EffectorP_secreted.gff/g')
    #ProgDir=/home/fanron/git_repos/tools/gene_prediction/ORF_finder
    #$ProgDir/extract_gff_for_sigP_hits.pl $OutFileHeaders $Gff effectorP ID > $EffectorP_Gff
  #done
  #Note: the Gff part is to know the location of effectors in genome.
```

V.alfafae - VaMs102
169

```bash
for File in $(ls analysis/effectorP/V.dahliae/JR2/*_EffectorP.txt); do
Test=$(echo "$File" | sed 's/_EffectorP.txt/_test.txt/'g)
cat "$File" | cut -d' ' -f1,7 | awk -v OFS="\t" '$1=$1' | cut -d$'\t' -f1,3,4 > $Test
Strain=$(echo $File | rev | cut -f2 -d '/' | rev)
Organism=$(echo $File | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
Headers=$(echo "$Test" | sed 's/_test.txt/_EffectorP_headers.txt/g')
cat $Test | grep 'Effector' | cut -f1 -d ' ' | awk -v OFS="\t" '$1=$1' | cut -d$'\t' -f1 > $Headers
Secretome=$(ls gene_pred/final_genes_signalp-4.1/V.dahliae/JR2/*_final_sp_no_trans_mem.aa)
OutFile=$(echo "$Test" | sed 's/_test.txt/_EffectorP_secreted.aa/g')
ProgDir=/home/fanron/git_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_from_fasta.py --fasta $Secretome --headers $Headers > $OutFile
OutFileHeaders=$(echo "$Test" | sed 's/_test.txt/_EffectorP_secreted_headers.txt/g')
cat $OutFile | grep '>' | tr -d '>' > $OutFileHeaders
cat $OutFileHeaders | wc -l
Gff=$(ls assembly/merged_canu_spades/V.dahliae/JR2/ensembl/Verticillium_dahliaejr2.GCA_000400815.2.33.gff3)
EffectorP_Gff=$(sed 's/_test.txt/_EffectorP_secreted.gff/g')
ProgDir=/home/fanron/git_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_gff_for_sigP_hits.pl $OutFileHeaders $Gff effectorP ID > $EffectorP_Gff
done
```
 | 
  V.dahliae - JR2
182

```bash
for File in $(ls analysis/effectorP/V.dahliae/Ls17/*_EffectorP.txt); do
Test=$(echo "$File" | sed 's/_EffectorP.txt/_test.txt/'g)
cat "$File" | cut -d' ' -f1,9 | awk -v OFS="\t" '$1=$1' | cut -d$'\t' -f1,3,4 > $Test
Strain=$(echo $File | rev | cut -f2 -d '/' | rev)
Organism=$(echo $File | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
Headers=$(sed 's/_test.txt/_EffectorP_headers.txt/g')
cat $File | grep 'Effector' | cut -f1 -d ' ' > $Headers
Secretome=$(ls gene_pred/final_genes_signalp-4.1/V.dahliae/Ls17/*_final_sp_no_trans_mem.aa)
OutFile=$(echo "$File" | sed 's/_test.txt/_EffectorP_secreted.aa/g')
ProgDir=/home/fanron/git_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_from_fasta.py --fasta $Secretome --headers $Headers > $OutFile
OutFileHeaders=$(echo "$File" | sed 's/_test.txt/_EffectorP_secreted_headers.txt/g')
cat $OutFile | grep '>' | tr -d '>' > $OutFileHeaders
cat $OutFileHeaders | wc -l
Gff=$(ls assembly/merged_canu_spades/V.dahliae/Ls17/ensembl/*.gff3)
EffectorP_Gff=$(sed 's/_test.txt/_EffectorP_secreted.gff/g')
ProgDir=/home/fanron/git_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_gff_for_sigP_hits.pl $OutFileHeaders $Gff effectorP ID > $EffectorP_Gff
done
```
  V.dahliae - Ls17
155

## C) CAZY proteins

Carbohydrte active enzymes were idnetified using CAZYfollowing recomendations
at http://csbl.bmb.uga.edu/dbCAN/download/readme.txt :

```bash
  for Proteome in $(ls assembly/merged_canu_spades/*/*/ensembl/*.pep.all.fa); do
    Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
    OutDir=gene_pred/CAZY/$Organism/$Strain
    mkdir -p $OutDir
    Prefix="$Strain"_CAZY
    CazyHmm=../../../dbCAN/dbCAN-fam-HMMs.txt
    ProgDir=/home/fanron/git_repos/tools/seq_tools/feature_annotation/HMMER
    # ProgDir=/home/gomeza/git_repos/emr_repos/tools/seq_tools/feature_annotation/HMMER
    qsub $ProgDir/sub_hmmscan.sh $CazyHmm $Proteome $Prefix $OutDir
  done
```
The Hmm parser was used to filter hits by an E-value of E1x10-5 or E 1x10-e3 if they had a hit over a length of X %.

********
Those proteins with a signal peptide were extracted from the list and gff files
representing these proteins made.

```bash
  for File in $(ls gene_pred/CAZY/V.dahliae/JR2/*CAZY.out.dm); do
      Strain=$(echo $File | rev | cut -f2 -d '/' | rev)
      Organism=$(echo $File | rev | cut -f3 -d '/' | rev)
      OutDir=$(dirname $File)
      echo "$Organism - $Strain"
      ProgDir=/home/groups/harrisonlab/dbCAN
      $ProgDir/hmmscan-parser.sh $OutDir/"$Strain"_CAZY.out.dm > $OutDir/"$Strain"_CAZY.out.dm.ps
      CazyHeaders=$(echo $File | sed 's/.out.dm/_headers.txt/g')
      cat $OutDir/"$Strain"_CAZY.out.dm.ps | cut -f3 | sort | uniq > $CazyHeaders
      echo "number of CAZY genes identified:"
      cat $CazyHeaders | wc -l
  done      
```
V.dahliae - JR2
number of CAZY genes identified:
606

```bash
for File in $(ls gene_pred/CAZY/V.dahliae/Ls17/*CAZY.out.dm); do
      Strain=$(echo $File | rev | cut -f2 -d '/' | rev)
      Organism=$(echo $File | rev | cut -f3 -d '/' | rev)
      OutDir=$(dirname $File)
      echo "$Organism - $Strain"
      ProgDir=/home/groups/harrisonlab/dbCAN
      $ProgDir/hmmscan-parser.sh $OutDir/"$Strain"_CAZY.out.dm > $OutDir/"$Strain"_CAZY.out.dm.ps
      CazyHeaders=$(echo $File | sed 's/.out.dm/_headers.txt/g')
      cat $OutDir/"$Strain"_CAZY.out.dm.ps | cut -f3 | sort | uniq > $CazyHeaders
      echo "number of CAZY genes identified:"
      cat $CazyHeaders | wc -l
      Gff=$(ls assembly/merged_canu_spades/V.*/Ls17/ensembl/*.gff3)
      CazyGff=$OutDir/"$Strain"_CAZY.gff
      ProgDir=/home/fanron/git_repos/tools/gene_prediction/ORF_finder
      $ProgDir/extract_gff_for_sigP_hits.pl $CazyHeaders $Gff CAZyme ID > $CazyGff

      SecretedProts=$(ls gene_pred/final_genes_signalp-4.1/$Organism/$Strain/"$Strain"_final_sp_no_trans_mem.aa)
      SecretedHeaders=$(echo $SecretedProts | sed 's/.aa/_headers.txt/g')
      cat $SecretedProts | grep '>' | tr -d '>' > $SecretedHeaders
      CazyGffSecreted=$OutDir/"$Strain"_CAZY_secreted.gff
      $ProgDir/extract_gff_for_sigP_hits.pl $SecretedHeaders $CazyGff Secreted_CAZyme ID > $CazyGffSecreted
      echo "number of Secreted CAZY genes identified:"
      cat $CazyGffSecreted | grep -w 'mRNA' | cut -f9 | tr -d 'ID=' | cut -f1 -d ';' > $OutDir/"$Strain"_CAZY_secreted_headers.txt
      cat $OutDir/"$Strain"_CAZY_secreted_headers.txt | wc -l
    done
```
V.dahliae - Ls17
number of CAZY genes identified:
623

```bash
for File in $(ls gene_pred/CAZY/*/VaMs102/*CAZY.out.dm); do
      Strain=$(echo $File | rev | cut -f2 -d '/' | rev)
      Organism=$(echo $File | rev | cut -f3 -d '/' | rev)
      OutDir=$(dirname $File)
      echo "$Organism - $Strain"
      ProgDir=/home/groups/harrisonlab/dbCAN
      $ProgDir/hmmscan-parser.sh $OutDir/"$Strain"_CAZY.out.dm > $OutDir/"$Strain"_CAZY.out.dm.ps
      CazyHeaders=$(echo $File | sed 's/.out.dm/_headers.txt/g')
      cat $OutDir/"$Strain"_CAZY.out.dm.ps | cut -f3 | sort | uniq > $CazyHeaders
      echo "number of CAZY genes identified:"
      cat $CazyHeaders | wc -l
      Gff=$(ls assembly/merged_canu_spades/V.*/VaMs102/ensembl/*.gff3)
      CazyGff=$OutDir/"$Strain"_CAZY.gff
      ProgDir=/home/fanron/git_repos/tools/gene_prediction/ORF_finder
      $ProgDir/extract_gff_for_sigP_hits.pl $CazyHeaders $Gff CAZyme ID > $CazyGff

      SecretedProts=$(ls gene_pred/final_genes_signalp-4.1/$Organism/$Strain/"$Strain"_final_sp_no_trans_mem.aa)
      SecretedHeaders=$(echo $SecretedProts | sed 's/.aa/_headers.txt/g')
      cat $SecretedProts | grep '>' | tr -d '>' > $SecretedHeaders
      CazyGffSecreted=$OutDir/"$Strain"_CAZY_secreted.gff
      $ProgDir/extract_gff_for_sigP_hits.pl $SecretedHeaders $CazyGff Secreted_CAZyme ID > $CazyGffSecreted
      echo "number of Secreted CAZY genes identified:"
      cat $CazyGffSecreted | grep -w 'mRNA' | cut -f9 | tr -d 'ID=' | cut -f1 -d ';' > $OutDir/"$Strain"_CAZY_secreted_headers.txt
      cat $OutDir/"$Strain"_CAZY_secreted_headers.txt | wc -l
    done
```
  V.alfafae - VaMs102
  number of CAZY genes identified:
  599

  ## D) Identify Small secreted cysteine rich proteins

  Small secreted cysteine rich proteins were identified within secretomes. These proteins may be identified by EffectorP, but this approach allows direct control over what constitutes a SSCP.
  ```bash
  for Secretome in $(ls gene_pred/final_genes_signalp-4.1/*/*/*_final_sp_no_trans_mem.aa); do
  Strain=$(echo $Secretome| rev | cut -f2 -d '/' | rev)
  Organism=$(echo $Secretome | rev | cut -f3 -d '/' | rev)
  echo "$Organism - $Strain"
  OutDir=analysis/sscp/$Organism/$Strain
  mkdir -p $OutDir
  ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/sscp
  $ProgDir/sscp_filter.py --inp_fasta $Secretome --max_length 400 --threshold 1 --out_fasta $OutDir/"$Strain"_sscp_all_results.fa
  cat $OutDir/"$Strain"_sscp_all_results.fa | grep 'Yes' > $OutDir/"$Strain"_sscp.fa
  echo "Number of effectors predicted by EffectorP:"
  EffectorP=$(ls analysis/effectorP/$Organism/$Strain/*_EffectorP_secreted_headers.txt)
  cat $EffectorP | wc -l
  echo "Number of SSCPs predicted by both effectorP and this approach"
  cat $OutDir/"$Strain"_sscp.fa | grep '>' | tr -d '>' > $OutDir/"$Strain"_sscp_headers.txt
  cat $OutDir/"$Strain"_sscp_headers.txt $EffectorP | cut -f1 | sort | uniq -d | wc -l
  done
  ```
V.alfafae - VaMs102
% cysteine content threshold set to:    1
maximum length set to:  300
No. short-cysteine rich proteins in input fasta:        266
Number of effectors predicted by EffectorP:
169
Number of SSCPs predicted by both effectorP and this approach
129
V.dahliae - 12008
% cysteine content threshold set to:    1
maximum length set to:  300
No. short-cysteine rich proteins in input fasta:        300
Number of effectors predicted by EffectorP:
187
Number of SSCPs predicted by both effectorP and this approach
149
V.dahliae - JR2
% cysteine content threshold set to:    1
maximum length set to:  300
No. short-cysteine rich proteins in input fasta:        284
Number of effectors predicted by EffectorP:
177
Number of SSCPs predicted by both effectorP and this approach
148
V.dahliae - Ls17
% cysteine content threshold set to:    1
maximum length set to:  300
No. short-cysteine rich proteins in input fasta:        272
Number of effectors predicted by EffectorP:
155
Number of SSCPs predicted by both effectorP and this approach
120































The number of SSCPs found in 12008 are: 4
The number of SSCPs found in JR2 are: 2
The number of SSCPs found in Ls17 are: 2
The number of SSCPs found in VaMs102 are: 1




















number of LysM and NPP1 in verticillium strains
cat ..interproscan..gff3 | grep 'LysM' | wc -l
................................'NPP1'........

     12008  JR2  Ls17  VaM102  51  53  58  61
LysM  11     11   11     11    11  11  11  11

NPP1   7     7    9      8     7   7   7    8
   