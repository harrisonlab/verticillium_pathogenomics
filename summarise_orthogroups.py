#!/usr/bin/python

'''
Summarise ortholog groups by
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
ap.add_argument('--set1',required=True,type=str,nargs='+',help='Set of isolates to compare vs group2')
ap.add_argument('--set2',required=True,type=str,nargs='+',help='Set of isolates to compare vs group1')
ap.add_argument('--all',required=True,type=str,nargs='+',help='All isolates used in orthology analysis')
ap.add_argument('--effectors',required=True,type=str,help='list of genes identified as Vd effectors')


conf = ap.parse_args()
grp1_isolates = conf.set1
grp2_isolates = conf.set2
all_isolates = conf.all

with open(conf.orthogroups) as f:
    orthogroup_lines = f.readlines()
with open(conf.effectors) as f:
    effector_lines = f.readlines()

effector_set = Set()
for line in effector_lines:
    line = line.rstrip("\n")
    effector_set.add(line)


# strain_id = strain_id + "|"

orthogroup_dict = defaultdict(list)
orthogroup_content_dict = defaultdict(list)
# grp1_isolates = ["Fus2", "125", "A23"]
# grp2_isolates = ["A28", "D2", "PG"]
# all_isolates = grp1_isolates + grp2_isolates

header_line = []
header_line.extend(["Orthogroup ID", "Presence status"])
header_line.extend(all_isolates)
header_line.extend(["Expansion", "Contained effectors", "Contained genes"])
print "\t".join(header_line)

for line in orthogroup_lines:
    line = line.rstrip("\n")
    split_line = line.split(" ")
    orthogroup_id = split_line[0].replace(":", "")
    orthogroup_contents = []
    orthogroup_content_dict.clear()
    for isolate in all_isolates:
        num_genes = line.count((isolate + "|"))
        # orthogroup_contents.append(str(isolate) + "(" + str(num_genes) + ")")
        orthogroup_contents.append(str(num_genes))
        content_str = "\t".join(orthogroup_contents)
        orthogroup_content_dict[isolate] = num_genes
    effector_list = []
    for gene in split_line[1:]:
        if gene in effector_set:
            effector_list.append(gene)
    grp1_numbers = []
    for isolate in grp1_isolates:
        grp1_numbers.append(orthogroup_content_dict[isolate])
    max_path = max(grp1_numbers)
    min_path = min(grp1_numbers)
    grp2_numbers = []
    for isolate in grp2_isolates:
        grp2_numbers.append(orthogroup_content_dict[isolate])
    max_non_path = max(grp2_numbers)
    min_non_path = min(grp2_numbers)
    if min_path > max_non_path:
        expansion_status = "grp1_expanded"
    elif min_non_path > max_path:
        expansion_status = "grp2_expanded"
    else:
        expansion_status = ""
    contained_effectors = ";".join(effector_list)
    contained_genes = ";".join(split_line[1:])
    if all(x in line for x in all_isolates):
        print("\t".join([orthogroup_id, "all_isolates", content_str, expansion_status, contained_effectors, contained_genes]))
    elif all(x not in line for x in grp2_isolates) and all(x in line for x in grp1_isolates):
        print("\t".join([orthogroup_id, "grp2_absent", content_str, expansion_status, contained_effectors, contained_genes]))
    elif all(x in line for x in grp2_isolates) and all(x not in line for x in grp1_isolates):
        print("\t".join([orthogroup_id, "grp1_absent", content_str, expansion_status, contained_effectors, contained_genes]))
    elif any(x in line for x in all_isolates):
        print("\t".join([orthogroup_id, "some_isolates", content_str, expansion_status, contained_effectors, contained_genes]))

    # for gene_id in split_line[1:]:
    #     if strain_id in gene_id:
    #         gene_id = gene_id.replace(strain_id, "")
    #         # orthogroup_dict[gene_id] = [orthogroup_id]
    #
    #         if all(x in line for x in all_isolates):
    #             orthogroup_dict[gene_id].extend([orthogroup_id, "all_isolates", content_str, expansion_status])
    #         elif all(x not in line for x in grp2_isolates) and all(x in line for x in grp1_isolates):
    #             orthogroup_dict[gene_id].extend([orthogroup_id, "grp1_isolates_all", content_str, expansion_status])
    #         elif all(x in line for x in grp2_isolates) and all(x not in line for x in grp1_isolates):
    #             orthogroup_dict[gene_id].extend([orthogroup_id, "grp2_isolates_all", content_str, expansion_status])
    #         elif any(x in line for x in all_isolates):
    #             orthogroup_dict[gene_id].extend([orthogroup_id, "some_isolates", content_str, expansion_status])
