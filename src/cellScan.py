#!/usr/bin/env python3

import argparse
import sys
import os
import logging
import shutil
import glob
import re
import pysam
import time
import subprocess
import pandas as pd
import numpy as np
import gzip
from pysam import VariantFile
import multiprocessing as mp
from multiprocessing import Pool



def robust_get_tag(read, tag_name):  
	try:  
		return read.get_tag(tag_name)
	except KeyError:
		return "NotFound"

def rev_compl(st):
    nn = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(nn[n] for n in reversed(st))

#para_lst = para.strip().split(":")
#chr = para_lst[0]
#cell = para_lst[1]
#out = para_lst[2]
#app_path = para_lst[3]

#assert os.path.isfile(bam_filter), "*.fiter.targeted.bam file {} cannot be found!".format(bam_filter)
in_fasta = "/rsrch3/scratch/bcb/jdou1/scAncestry/ref/fasta/genome.fa"
#in_vcf = "/rsrch3/scratch/bcb/jdou1/bam_process/chr1.subset.vcf"
#in_bam = "/rsrch3/scratch/bcb/jdou1/bam_process/chr1.filter.targeted.1M.bam"
in_vcf = "/rsrch3/scratch/bcb/jdou1/bam_process/chr1.gp.vcf"
in_bam = "/rsrch3/scratch/bcb/jdou1/bam_process/chr1.filter.targeted.bam"


ref_fa = pysam.FastaFile(in_fasta)
vcf = VariantFile(in_vcf) 
snv_info = {}
index = 0 
motif_len = 3

for rec in vcf.fetch():
	seq = ref_fa.fetch(rec.chrom, rec.pos-4, rec.pos+3)
	seq_rev_compl = rev_compl(seq)
	id = str(rec.chrom)+":"+str(rec.pos) + ":" + rec.ref + ":" + rec.alts[0]
	mydict = dict(chr = rec.chrom, pos = rec.pos, motif_pos= seq, ref_allele = rec.ref, alt_allele=rec.alts[0])
	#print(id + " " + seq + " " + seq_rev_compl)

	index = index + 1
	snv_info[index]=mydict  

snv_tol = index 
index = 1
lower_index = 1


read_tol = 0 
read_cover = 0
read_wild = 0 
read_mutated = 0
read_noAllele = 0
infile = pysam.AlignmentFile(in_bam,"rb")
fp = open("search.log", "wt")
for s in infile:
	t  = robust_get_tag(s,"CB")
	read_name = s.query_name 
	align_chr = s.reference_name
	align_start = s.reference_start
	align_seq = s.query_sequence

	mystart = str(snv_info[lower_index]["pos"])
	read_tol = read_tol + 1 
	read_len = len(align_seq)
	# print("--------> read  started:" + str(align_start)  + "  snv:" + str(snv_info[index_start]["pos"]))

	### get the SNVs locating in [align_start, align_start + read_len] ### 
	if read_tol%1000000 == 0:
		print("scanning read " + str(read_tol)) 

	lock = 0 
	flag = "move"
	snv_cover = ""
	for i in range(lower_index,snv_tol,1):
		snv_pos = snv_info[i]["pos"]
		if snv_pos >= align_start:
			if lock==0: 
				lower_index = i 
			lock = lock + 1 
			if snv_pos <= align_start + read_len: 
				read_cover = read_cover + 1 
				snv_cover = str(snv_cover) + ":" + str(snv_pos)
				motif_pos = snv_info[i]["motif_pos"]
				motif_neg = motif_pos[:motif_len] + snv_info[i]["alt_allele"] + motif_pos[motif_len + 1:]
			
				if re.search(motif_pos, align_seq):
					read_wild = read_wild + 1 
				elif re.search(motif_neg, align_seq):
					read_mutated = read_mutated + 1
				else:
					delta = snv_pos - align_start	
					if read_len - delta < motif_len:
						motif_pos_part = motif_pos[0:motif_len+1]
						motif_neg_part = motif_neg[0:motif_len+1]
						seq_part = align_seq[read_len-2*motif_len-1:read_len]
						
						if re.search(motif_pos_part, seq_part):
							read_wild = read_wild + 1 			
						elif re.search(motif_neg_part, seq_part):	
							read_mutated = read_mutated + 1

					elif delta <= motif_len: 
						motif_pos_part = motif_pos[motif_len:len(motif_pos)]
						motif_neg_part = motif_neg[motif_len:len(motif_neg)]
						seq_part = align_seq[0:2*motif_len+1]

						if re.search(motif_pos_part, seq_part):
							read_wild = read_wild + 1 
						elif re.search(motif_neg_part, seq_part):
							read_mutated = read_mutated + 1

					else:
						fp.write(str(snv_info[i])+"\n")
						fp.write(str(align_chr) + ":" + str(align_start) + ":" + align_seq+"\n") 
						read_noAllele = read_noAllele + 1 
			else:
				break 
        
infile.close()
fp.close()
print(str(read_tol) + ":" + str(read_cover) + ":" + str(read_wild) + ":" + str(read_mutated) + ":" + str(read_noAllele))
