#!/usr/bin/env python3
"""
The main interface of scPopGene
"""

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
from pysam import VariantFile

LIB_PATH = os.path.abspath(
	os.path.join(os.path.dirname(os.path.realpath(__file__)), "pipelines/lib"))

if LIB_PATH not in sys.path:
	sys.path.insert(0, LIB_PATH)

PIPELINE_BASEDIR = os.path.dirname(os.path.realpath(sys.argv[0]))
CFG_DIR = os.path.join(PIPELINE_BASEDIR, "cfg")

#import pipelines
#from pipelines import get_cluster_cfgfile
#from pipelines import PipelineHandler

# global logger
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter(
	'[{asctime}] {levelname:8s} {filename} {message}', style='{'))
logger.addHandler(handler)



def addChr(args):
	# edit the sequence names for your output header
	in_bam = args.bamFile
	prefix = 'chr'
	out_bam=in_bam+"tmp.bam"
	print(in_bam)
	input_bam = pysam.AlignmentFile(in_bam,"rb")
	new_head = input_bam.header.to_dict()
	for seq in new_head['SQ']:
		seq['SN'] = prefix  + seq['SN']
	# create output BAM with newly defined header
	with pysam.AlignmentFile(out_bam, "wb", header=new_head) as outf:
		for read in input_bam.fetch():
			prefixed_chrom = prefix + read.reference_name
			a = pysam.AlignedSegment(outf.header)
			a.query_name = read.query_name
			a.query_sequence = read.query_sequence
			a.reference_name = prefixed_chrom
			a.flag = read.flag
			a.reference_start = read.reference_start
			a.mapping_quality = read.mapping_quality
			a.cigar = read.cigar
			a.next_reference_id = read.next_reference_id
			a.next_reference_start = read.next_reference_start
			a.template_length = read.template_length
			a.query_qualities = read.query_qualities
			a.tags = read.tags
			outf.write(a)
	input_bam.close()
	outf.close()
	os.system("samtools  index " +  out_bam)
	os.system(" mv " + out_bam + " " + in_bam)
	os.system(" mv " + out_bam + ".bai  " + in_bam + ".bai")




def sort_chr(chr_lst):
	# sort chr IDs from 1...22
	chr_lst_sort = []
	for i in range(1, 23):
		i = str(i)
		if  i in chr_lst:
			chr_lst_sort.append(i)
		i_chr = "chr"+i 
		if  i_chr in chr_lst:
			chr_lst_sort.append(i_chr)
	chr_lst = chr_lst_sort 
	return chr_lst




def bamSplit(para):
	para_lst = para.strip().split(":")
	chr = para_lst[0]
	cell = para_lst[1]
	out = para_lst[2]
	app_path = para_lst[3]

	#assert os.path.isfile(bam_filter), "*.fiter.targeted.bam file {} cannot be found!".format(bam_filter)
	samtools = app_path + "/samtools" 
	output_bam =  out + "/Bam/merge.filter.targeted.bam" 
	infile = pysam.AlignmentFile(output_bam,"rb")
	# Note to change the read groups 
	tp =infile.header.to_dict()
	if len(tp['RG'])>1:
		tp['RG']= [tp['RG'][0]]
	tp['RG'][0]['SM'] = cell
	tp['RG'][0]['ID'] = cell
	cnt = 0
	outfile =  pysam.AlignmentFile( out + "/Bam/split_bam/" + cell  + ".bam", "wb", header=tp)
	for s in infile:
		t  = robust_get_tag(s,"CB")
		if not t=="NotFound":
			if t==cell:
				outfile.write(s)
				cnt = cnt +1 
		else:
			if re.search(cell, s.query_name):
				outfile.write(s)
				cnt = cnt + 1
	
	outfile.close()
	infile.close()
	cmd=samtools + " index " + out + "/Bam/split_bam/" + cell + ".bam"

	os.system(samtools + " index " + out + "/Bam/split_bam/" + cell + ".bam")

	return(cnt)




def jointCall(para):

	para_lst = para.strip().split(">")
	jobid = para_lst[0]
	chr=para_lst[1]
	out = para_lst[2]
	app_path = para_lst[3]
	reference = para_lst[4]
	samtools = app_path + "/samtools" 
	bcftools = app_path + "/bcftools" 
	bgzip = app_path + "/bgzip"
	bam_filter = out + "/Bam/split_bam/cell.bam.lst"
	cmd1 = samtools + " mpileup -b " + bam_filter + " -f "  + reference + " -r " +  jobid + " -q 20 -Q 20 -t DP4 -d 10000 -v "
	cmd1 = cmd1 + " | " + bcftools + " view " + " | "  + bcftools  + " norm -m-both -f " + reference 
	cmd1 = cmd1 + " | grep -v \"<X>\" | grep -v INDEL |" + bgzip +   " -c > " + out + "/somatic/" +  jobid + ".cell.gl.vcf.gz" 		
	print(cmd1)
	output = os.system(cmd1)

	# delete the bam files once snv calling was finished in specfic regions
	f = open(bam_filter, "r")
	for x in f:
		x = x.strip()
		#os.system("rm " + x)
		#os.system("rm " + x + ".bai")		
	f.close()

	if output == 0:
		return(jobid)
