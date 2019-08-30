#!/usr/bin/python

'''
UTR regions are publicaly available for the JR2 genome.
This script extracts -1000/+50 of the start of 5'-UTR or the 1500 bp
upstream of a gene if a UTR feature is missing.
If another gene is within 1000 bp then a shortened sequence
is predicted.
'''

# for feature in gff if gene then extract start/stop
# if followed by UTR region then extract start/stop
# make a dictionary of gene boundaries and contig boundaries.
# for each gene extract the sequence -1000/+50 of TSS
#     or -1000+0 of gene start site
# If another gene is within 1000bp identify this by selecting a range of
# 1000 values less than the TSS and see if they are in the dictionary of
# gene boundaries (reverse the list).
# Adjust to the number in the dictionary +1.
# reverse complement the sequence if gene is on the reverse strand.
# Save sequence to list
#
# print list.

import sys,argparse
import re
from collections import defaultdict
from sets import Set
from Bio.Seq import Seq
from Bio import SeqIO
import numpy as np

#-----------------------------------------------------
# Step 1
# Import variables & load input files
#-----------------------------------------------------
ap = argparse.ArgumentParser()
ap.add_argument('--gff',required=True,type=str,help='input gff file')
ap.add_argument('--fasta', required = False, type=str, default = False, help = 'File name to store converted gene names if required')
ap.add_argument('--prefix', required = False, type=str, default = False, help = 'outfile prefix')
ap.add_argument('--distance', required = False, type=int, default = False, help = 'distance extracted upstream of TSS (+500 if no UTR predicted5)')
conf = ap.parse_args() #sys.argv

# with open(conf.fasta) as f:
#     fasta_lines = f.readlines()
# for rec in SeqIO.parse(conf.fasta,"fasta"):
# 	seq = str(rec.seq)
record_dict = SeqIO.to_dict(SeqIO.parse(conf.fasta, "fasta"))

with open(conf.gff) as f:
    gff_lines = f.readlines()

prefix = conf.prefix
dist = conf.distance

#-----------------------------------------------------
# Step 2
# Define classes
#-----------------------------------------------------

class Feature_obj(object):
    def __init__(self):
        """Return a Annot_obj whose name is *transcript_id*"""
        self.id = ''
        self.feature_dict = defaultdict(list)
        self.TSS_dist = ''
        # self.start = ''
        # self.stop = ''
        # self.strand = ''
    def set_TSS_distance(self):
        if self.feature_dict["gene"]:
            TSS_obj = self.feature_dict["gene"][0]
            gene_obj = self.feature_dict["CDS"][0]
            if gene_obj.strand == '+':
                self.TSS_dist = int(gene_obj.start) - int(TSS_obj.start)
            elif gene_obj.strand == '-':
                self.TSS_dist = int(TSS_obj.stop) - int(gene_obj.stop)


class Gff_obj(object):
    # """A gene identified as differentially expressed in one isolate.
    # Attributes:
    #     transcript_id: A string representing the gene name.
    #     conditions_tested: List of conditions that the gene was tested for DE.
    #     conditions_positive: A binary list of integers representing whether the
    #         gene tested positive for each condition listed in conditions tested.
    # """

    def __init__(self):
        """Return a Annot_obj whose name is *transcript_id*"""
        self.contig = ''
        self.source = ''
        self.feature = ''
        self.start = ''
        self.stop = ''
        self.strand = ''
        self.evidence = ''
        self.codon = ''
        # self.attributes = ''
        self.id = ''
        self.parent = ''
        self.name = ''

    def set_conditions(self, gff_line):
        """Reset conditions_tested to a list of conditions"""
        gff_elements = gff_line.split("\t")
        self.contig = gff_elements[0]
        self.source = gff_elements[1]
        self.feature = gff_elements[2]
        self.start = gff_elements[3]
        self.stop = gff_elements[4]
        self.evidence = gff_elements[5]
        self.strand = gff_elements[6]
        self.codon = gff_elements[7]
        # col9_list = gff_elements[8].split(";")
        # print(col9_list)
        id = re.findall(r"Dbxref=GeneID:.+?[,;]", gff_elements[8])
        parent = re.findall(r"Parent=.*", gff_elements[8])
        name = re.findall(r"Name=.+?;", gff_elements[8])
        self.id = ",".join(id).replace("Dbxref=", "").replace(",", "").replace(";", "").replace(":", "_").replace(' ', "_")
        self.parent = ",".join(parent).replace(";", "").replace("Parent=", "").replace(' ', "_")
        self.name = ",".join(name).replace("Name=", "").replace(";", "").replace(' ', "_")



#-----------------------------------------------------
# Step 3
# Read GFF
#-----------------------------------------------------

prev_id = ""
obj_dict = defaultdict(list)
for line in gff_lines:
    if line.startswith('#'):
        continue
    line = line.rstrip()
    # print line
    split_line = line.split("\t")
    feature = split_line[2]
    if feature == 'tRNA':
        feature = 'CDS'
    if not any(feature == x for x in ['gene', 'CDS', 'five_prime_UTR']):
        continue
    col_9 = split_line[8]
    id = re.findall(r"Dbxref=GeneID:.+?[,;]", col_9)
    # print(id)
    # exit()
    if id:
        this_gene = id[0].replace("Dbxref=", "").replace(",", "").replace(";", "").replace(":", "_").replace(' ', "_")
    if not id:
        id = prev_id
    # this_gene = re.sub(r"\..*","", this_gene)
    prev_id = id
    # print("\t".join([this_gene, feature]))

    # Prepare a feature object:
    feat_obj = Gff_obj()
    feat_obj.set_conditions(line)
    # print feat_obj.contig
    # print "\t".join([feat_obj.contig, feat_obj.feature, feat_obj.id, feat_obj.parent])

    # Add the feature obj to the object dictionary:
    if obj_dict[this_gene]:
        obj_dict[this_gene].feature_dict[feature].append(feat_obj)
    else:
        gene_obj = Feature_obj()
        obj_dict[this_gene] = gene_obj
        obj_dict[this_gene].feature_dict[feature].append(feat_obj)

# print this_gene
gene_mRNA_objs = obj_dict[this_gene].feature_dict['CDS']
gene_mRNA_obj = gene_mRNA_objs[0]
# print "\t".join([gene_mRNA_obj.contig, gene_mRNA_obj.feature, gene_mRNA_obj.id, gene_mRNA_obj.parent])



#-----------------------------------------------------
# Step 4
# Create a dictionary of gene boundaries
#-----------------------------------------------------

gene_boundary_dict = defaultdict(set)
for key in obj_dict.keys():
    gene_obj = obj_dict[key]
    if gene_obj.feature_dict["CDS"]:
        for feat_obj in gene_obj.feature_dict["CDS"]:
            gene_boundary_dict[feat_obj.contig].add(feat_obj.start)
            gene_boundary_dict[feat_obj.contig].add(feat_obj.stop)
    # if gene_obj.feature_dict["CDS"]:
    #     feat_obj = gene_obj.feature_dict["CDS"][0]
    #     gene_boundary_dict[feat_obj.contig].add(feat_obj.start)
    #     gene_boundary_dict[feat_obj.contig].add(feat_obj.stop)

# print(feat_obj.contig)
# print(len(gene_boundary_dict[feat_obj.contig]))

for key in record_dict.keys():
    gene_boundary_dict[key].add("1")
    gene_boundary_dict[key].add("0")
    # print(len(record_dict[key]))
    gene_boundary_dict[key].add(len(record_dict[key]))
    gene_boundary_dict[key].add(len(record_dict[key])+1)

# print(record_dict[key].seq[1:100])

#-----------------------------------------------------
# Step 5
# For each gene extract the promoter region
#-----------------------------------------------------

out_feat_dict = defaultdict(list)
gff_outlines = ['#gff3']
fasta_outlines = []
for key in obj_dict.keys():
    gene_obj = obj_dict[key]
    if gene_obj.feature_dict['gene'] and gene_obj.feature_dict['CDS']:
        feat_obj = gene_obj.feature_dict["gene"][0]
        strand = feat_obj.strand
        start = int(feat_obj.start)
        stop = int(feat_obj.stop)
        contig = feat_obj.contig
        # TSS not present for gene
        # print key
        # print(gene_obj.feature_dict.keys())
        # if gene_obj.feature_dict['five_prime_UTR']:
            # print gene_obj.feature_dict['five_prime_UTR']
            # print "monkeys"
        feat_obj = gene_obj.feature_dict["gene"][0]
        strand = feat_obj.strand
        start = int(feat_obj.start)
        stop = int(feat_obj.stop)
        contig = feat_obj.contig
        name = gene_obj.feature_dict["gene"][0].name
        # print name
        if feat_obj.strand == "+":
            start = int(gene_obj.feature_dict["CDS"][0].start)
            promoter_start = int(start)
            for i in reversed(range(start - (dist), start, 1)):
                if str(i) in gene_boundary_dict[contig]:
                    # print("monkeys")
                    promoter_start = i
                    break
                promoter_start = i
            promoter_stop = start -1
            seq = record_dict[contig].seq[promoter_start:promoter_stop]
        elif feat_obj.strand == "-":
            # print gene_obj.feature_dict["gene"][0].id
            stop = int(gene_obj.feature_dict["CDS"][0].stop)
            promoter_stop = int(stop)
            for i in range(stop + 1, stop + (1 + dist), 1):
                if str(i) in gene_boundary_dict[contig]:
                    # print("monkeys")
                    promoter_stop = i
                    break
                promoter_stop = i
            promoter_start = stop + 1
            seq = record_dict[contig].seq[promoter_start:promoter_stop]
            seq.reverse_complement()
        out_key = "_".join([contig, str(promoter_start), str(promoter_stop)])
        if out_feat_dict[out_key]:
            prev_name = out_feat_dict[out_key][0].split("=")[-1]
            out_feat_dict[out_key][0] = "\t".join([contig, "extract_promoters", "IG_region", str(promoter_start), str(promoter_stop), '.', '.', '.', 'ID=' + prev_name + "_" + name])
            if len(seq) != 0:
                out_feat_dict[out_key][1] = "\n".join([">" + prev_name + "_" + name + "_IG_region-" + str(dist), str(seq)])
        else:
            out_feat_dict[out_key] = ["\t".join([contig, "extract_promoters", "IG_region", str(promoter_start), str(promoter_stop), '.', strand, '.', 'ID=' + name])]
            if len(seq) != 0:
                out_feat_dict[out_key].append("\n".join([">" + name + "_IG_region-" + str(dist), str(seq)]))
# print "\n".join(gff_outlines)

for key in out_feat_dict.keys():
    gff = out_feat_dict[key][0]
    # print out_feat_dict[key]
    if len(out_feat_dict[key]) == 2:
        seq = out_feat_dict[key][1]
    else:
        # print "monkeys"
        print gff
    gff_outlines.append(gff)
    fasta_outlines.append(seq)

outgff = ".".join([prefix, "gff"])
f = open(outgff,"w")
f.write("\n".join(gff_outlines))
f.close()

# print "\n".join(fasta_outlines)
outfasta = ".".join([prefix, "fa"])
f = open(outfasta,"w")
f.write("\n".join(fasta_outlines))
f.close()


# Retreive TSS distances from start codons
# TSS_distances = []
# for key in obj_dict.keys():
#     obj_dict[key].set_TSS_distance()
#     dist = obj_dict[key].TSS_dist
#     if dist != '':
#         TSS_distances.append(dist)

# print("Average distance of TSS from start codon:")
# print(np.mean(TSS_distances))


# outdist = "_".join([prefix, "TSS_distances.txt"])
# f = open(outdist,"w")
# f.write("\n".join(str(x) for x in TSS_distances))
# f.close()



# print contig
# print sorted(gene_boundary_dict[contig])
        # gene_boundary_dict[feat_obj.contig].add(feat_obj.start)
        # gene_boundary_dict[feat_obj.contig].add(feat_obj.stop)
