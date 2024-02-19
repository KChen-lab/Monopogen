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

# Library paths
LIB_PATH = os.path.abspath(
	os.path.join(os.path.dirname(os.path.realpath(__file__)), "pipelines/lib"))

if LIB_PATH not in sys.path:
	sys.path.insert(0, LIB_PATH)

# Pipeline paths and configurations - is this necessary?
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

# function to provide parameters given
def print_parameters_given(args):
	logger.info("Parameters in effect:")
	for arg in vars(args):
		if arg=="func": continue
		logger.info("--{} = [{}]".format(arg, vars(args)[arg]))

# function to check the existence of sample list file
def validate_sample_list_file(args):
	if args.check_hard_clipped:
		out=os.popen("command -v bioawk").read().strip()
		assert out!="", "Program bioawk cannot be found!"

	assert os.path.isfile(args.sample_list), "Sample index file {} cannot be found!".format(args.sample_list)

	try:
		with open(args.sample_list) as f_in:
			for line in f_in:
				record = line.strip().split("\t")
				logger.debug("Checking sample {}".format(record[0]))
				assert len(record)==3, "Every line has to have exactly 3 tab-delimited columns! Line with sample name {} does not satisify this requiremnt!".format(record[0])
				assert os.path.isfile(record[1]), "Bam file {} cannot be found!".format(record[1])
				assert os.path.isfile(record[1]+".bai"), "Bam file {} has not been indexed!".format(record[1])
				assert os.path.isabs(record[1]), "Please use absolute path for bam file {}!".format(record[1])

				if args.check_hard_clipped:
					logger.debug("Checking existence of hard-clipped reads.")
					cmd = "samtools view {} | bioawk -c sam 'BEGIN {{count=0}} ($cigar ~ /H/)&&(!and($flag,256)) {{count++}} END {{print count}}'".format(record[1])
					logger.debug("Command: "+cmd)
					out=os.popen(cmd).read().strip()
					logger.debug("Results: "+out)
					assert out=="0", "Bam file {} contains hard-clipped reads without proper flag (0x100) set! Please use -M or -Y options of BWA MEM!".format(record[1])

				try:
					float(record[2])
					assert 0.0 <= float(record[2]) and float(record[2]) <= 1.0, "Contamination rate of sample {0} has to be a float number between 0 and 1 instead of {1}!".format(record[0], record[2])
				except:
					logger.error("Contamination rate of sample {0} has to be a float number between 0 and 1 instead of {1}!".format(record[0], record[2]))
					exit(1)

	except Exception:
		logger.error("There is something wrong with the sample index file. Check the logs for more information.")
		print(sys.exc_info())
            
		raise sys.exc_info()[0]

# function to check the existence of reference genome
def validate_user_setting_germline(args):
	#assert os.path.isfile(args.bamFile), "The bam list file {} cannot be found!".format(args.bamFile)
	assert os.path.isfile(args.reference), "The genome reference fasta file {} cannot be found!".format(args.reference)
	assert os.path.isdir(args.imputation_panel), "Filtered genotype file of 1KG3 ref panel {} cannot be found!".format(args.imputation_panel)
	assert os.path.isfile(args.region), "The region file {} cannot be found!".format(args.region)
	# check whether each bam file available	
	for chr in range(1, 23):
		bamFile = args.out + "/Bam/chr" +  str(chr) +  ".filter.bam.lst"
		with open(bamFile) as f_in:
			for line in f_in:
				line = line.strip()
				assert os.path.isfile(line), "The bam file {} cannot be found!".format(line)
				assert os.path.isfile(line + ".bai"), "The bam.bai file {} cannot be found!".format(line)
	# check whether region files were set correctly 
	with open(args.region) as f_in:
		for line in f_in:
			record = line.strip().split(",")
			assert len(record)==3 or len(record)==1, "Every line has to have exactly 3 comma-delimited columns chr1,1,100000 or chr1 (on the whole chromosome)! Line with region {} does not satisify this requiremnt!".format(line)

# function to check the existence of tools
def check_dependencies(args):
	# check whether each program is available
    programs_to_check = ["beagle.27Jul16.86a.jar"]
    system_programs_to_check = ["vcftools", "bcftools", "samtools", "bgzip", "java"]

    print(f"Check existence of the proper tools ({programs_to_check} & {system_programs_to_check})...")
    # check whether each program is available
    for prog in programs_to_check:
        if args.verbose:
            print(f"> checking existence of [{prog}]")
        prog_path = os.path.join(args.app_path, prog)
        assert os.path.exists(prog_path), f"Program {prog} cannot be found at {prog_path}!"

    for prog_system in system_programs_to_check:
        if args.verbose:
            print(f"> checking existence of [{prog_system}]")
        out_system = os.popen(f"command -v {prog_system}").read()
        assert out_system.strip() != "", f"Program {prog_system} cannot be found!"

# function to align the bam file
def addChr(in_bam, samtools):
	# edit the sequence names for your output header
	prefix = 'chr'
	out_bam=in_bam+"tmp.bam"
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
	os.system(samtools + " index " +  out_bam)
	os.system(" mv -v " + out_bam + " " + in_bam)
	os.system(" mv -v " + out_bam + ".bai  " + in_bam + ".bai")

# function to filter the bam file
def BamFilter(myargs):
	bamFile = myargs.get("bamFile")
	search_chr = myargs.get("chr")
	samtools = myargs.get("samtools")
	chr = search_chr
	id = myargs.get("id")
	#samtools = myargs.get("samtools")
	max_mismatch = myargs.get("max_mismatch")
	out = myargs.get("out")

	# check whether the output directory exists - is this necessary, as it is already checked in the Monopogen.py?
	os.system("mkdir -p -v " + out +  "/Bam")
	infile = pysam.AlignmentFile(bamFile,"rb")
	contig_names = infile.references
	cnt=0 
	for contig in contig_names:
		if contig.startswith("chr"):
			cnt=cnt+1
	if cnt==0:
			logger.info("The contig {} does not contain the prefix 'chr' and we will add 'chr' on it ".format(search_chr))
			search_chr = search_chr[3:]
			#searchBam = args.bamFile+".bam"
			#addChr(args.bamFile, searchBam)
			#infile.close()
			#infile = pysam.AlignmentFile(searchBam,"rb")

	tp =infile.header.to_dict()
			
	#print(tp)
	#To avoid the format issue, we update the RG flag based on sample information

	#if not "RG" in tp:
	sampleID = os.path.splitext(os.path.basename(myargs["bamFile"]))[0]
	tp1 = [{'SM':sampleID,'ID':sampleID, 'LB':"0.1", 'PL':"ILLUMINA", 'PU':sampleID}]
	tp.update({'RG': tp1})
		#print(tp)
		#tp['RG'][0]['SM'] = 
		#tp['RG'][0]['ID'] = os.path.splitext(os.path.basename(args.bamFile))[0]
  
	outfile =  pysam.AlignmentFile( out + "/Bam/" +id + "_" + chr + ".filter.bam", "wb", header=tp)


	for s in infile.fetch(search_chr):  
		#print(str(s.query_length)  + ":" + str(s.get_tag("AS")) + ":" + str(s.get_tag("NM")))
		if s.has_tag("NM"):
			val= s.get_tag("NM")
		if s.has_tag("nM"):
			val= s.get_tag("nM")                  
		if val < max_mismatch:
			outfile.write(s)
	infile.close()
	outfile.close()

	os.system(samtools + " index " +  out + "/Bam/" + id+ "_"  + chr + ".filter.bam")
	if cnt ==0:
		addChr(out + "/Bam/" +  id+ "_" + chr+ ".filter.bam", samtools)
	bamfile = out + "/Bam/" +  id+ "_" + chr+ ".filter.bam"
	return(bamfile)
	#args.bam_filter = args.out + "/Bam/" + args.chr + ".filter.bam"


def robust_get_tag(read, tag_name):  
	try:  
		return read.get_tag(tag_name)
	except KeyError:
		return "NotFound"

def runCMD(cmd):
	output = os.system(cmd)
	if output == 0:
		return(region)
	#process = subprocess.run(cmd, shell=True, stdout=open(args.logfile, 'w'), stderr=open(args.logfile,'w'))
