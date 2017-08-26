#!/usr/bin/python

'''
Orthology analysis was performed on Vd gene models prior to being submitted to ncbi
and subsequently renamed. This script renames the genes present in each
orthogroup based upon a conversion table generated as part of the genome submission
script.
'''


#-----------------------------------------------------
# Step 1
# Import variables & load input files
#-----------------------------------------------------

import sys
import argparse
import re
from sets import Set
from collections import defaultdict
from operator import itemgetter

ap = argparse.ArgumentParser()

ap.add_argument('--orthogroups',required=True,type=str,help='A fasta file of the assembled contigs')
ap.add_argument('--tsv_12008',required=True,type=str,help='A table of gene conversions')
ap.add_argument('--tsv_51',required=True,type=str,help='A table of gene conversions')
ap.add_argument('--tsv_53',required=True,type=str,help='A table of gene conversions')
ap.add_argument('--tsv_58',required=True,type=str,help='A table of gene conversions')
ap.add_argument('--tsv_61',required=True,type=str,help='A table of gene conversions')

conf = ap.parse_args()

with open(conf.orthogroups) as f:
    orthogroup_lines = f.readlines()

with open(conf.tsv_12008) as f:
    tsv_12008 = f.readlines()
with open(conf.tsv_51) as f:
    tsv_51 = f.readlines()
with open(conf.tsv_53) as f:
    tsv_53 = f.readlines()
with open(conf.tsv_58) as f:
    tsv_58 = f.readlines()
with open(conf.tsv_61) as f:
    tsv_61 = f.readlines()

#-----------------------------------------------------
# Step 2
# build dictionaries for renaming genes
#-----------------------------------------------------

conversion_dict = defaultdict(list)
for line in tsv_12008:
    line = line.rstrip()
    split_line = line.split()
    old_name = split_line[0]
    old_name = re.sub('^.+?_', '', old_name, 1) # replaces first occurence
    new_name = split_line[1]
    new_name = re.sub('^.+?_', '', new_name, 1)
    new_name = re.sub('^vAg', 'g', new_name)
    conversion_dict["12008|" + old_name] = "12008|" + new_name

for line in tsv_51:
    line = line.rstrip()
    split_line = line.split()
    old_name = split_line[0]
    old_name = re.sub('^.+?_', '', old_name, 1)
    new_name = split_line[1]
    new_name = re.sub('^.+?_', '', new_name, 1)
    conversion_dict["51|" + old_name] = "12251|" + new_name

for line in tsv_53:
    line = line.rstrip()
    split_line = line.split()
    old_name = split_line[0]
    old_name = re.sub('^.+?_', '', old_name, 1)
    new_name = split_line[1]
    new_name = re.sub('^.+?_', '', new_name, 1)
    conversion_dict["53|" + old_name] = "12253|" + new_name

for line in tsv_58:
    line = line.rstrip()
    split_line = line.split()
    old_name = split_line[0]
    old_name = re.sub('^.+?_', '', old_name, 1)
    new_name = split_line[1]
    new_name = re.sub('^.+?_', '', new_name, 1)
    conversion_dict["58|" + old_name] = "12158|" + new_name

for line in tsv_61:
    line = line.rstrip()
    split_line = line.split()
    old_name = split_line[0]
    old_name = re.sub('^.+?_', '', old_name, 1)
    new_name = split_line[1]
    new_name = re.sub('^.+?_', '', new_name, 1)
    conversion_dict["61|" + old_name] = "12161|" + new_name

for line in orthogroup_lines:
    line = line.rstrip()
    split_line = line.split()
    parsed_line = []
    for gene in split_line:
        split_gene = gene.split('.') # transcript id not in dictionary g1.t1
        if conversion_dict[split_gene[0]]:
            split_gene[0] = conversion_dict[split_gene[0]]
            gene = ".".join(split_gene)
            # print gene
        parsed_line.append(gene)
    print "\t".join(parsed_line)
