#!/usr/bin/env python3
"""
The main interface of monopgen
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
import numpy as np
import gzip
from pysam import VariantFile
import multiprocessing as mp
from multiprocessing import Pool


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




def print_parameters_given(args):
	logger.info("Parameters in effect:")
	for arg in vars(args):
		if arg=="func": continue
		logger.info("--{} = [{}]".format(arg, vars(args)[arg]))


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


def check_dependencies(args):
	programs_to_check = ("vcftools", "bgzip",  "bcftools", "beagle.08Feb22.fa4.jar", "beagle.27Jul16.86a.jar","samtools","picard.jar", "java")

	for prog in programs_to_check:
		out = os.popen("command -v {}".format(args.app_path + "/" + prog)).read()
		assert out != "", "Program {} cannot be found!".format(prog)

#	python_pkgs_to_check = ("drmaa",)

#	for pkg in python_pkgs_to_check:
#		out_pipe = os.popen('python -c "import {}"'.format(pkg))

#		assert out_pipe.close() is None, "Python module {} has not been installed!".format(pkg)



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
	os.system(" mv " + out_bam + " " + in_bam)
	os.system(" mv " + out_bam + ".bai  " + in_bam + ".bai")



def BamFilter(myargs):
	bamFile = myargs.get("bamFile")
	search_chr = myargs.get("chr")
	samtools = myargs.get("samtools")
	chr = search_chr
	id = myargs.get("id")
	#samtools = myargs.get("samtools")
	max_mismatch = myargs.get("max_mismatch")
	out = myargs.get("out")

	os.system("mkdir -p " + out +  "/Bam")
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

	if not "RG" in tp:
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
		addChr(out + "/Bam/" +  id+ "_" + chr+ ".filter.bam")
	bamfile = out + "/Bam/" +  id+ "_" + chr+ ".filter.bam"
	return(bamfile)
	#args.bam_filter = args.out + "/Bam/" + args.chr + ".filter.bam"

def getDPinfo(args):
	out = args.out
	gp_vcf_in = VariantFile(out + "/germline/" +  args.chr + ".gp.vcf.gz") 
	gp_info = {}
	gp_info_GT = {}
	error_rate = {}
	error_rate_cnt = {}
	for rec in gp_vcf_in.fetch():
		GT = [value['GT'] for value in rec.samples.values()][0]
		GP = [value['GP'] for value in rec.samples.values()][0]
		GT = str(GT[0]) + "/" + str(GT[1])
		GP = str(GP[0]) + "," + str(GP[1]) + "," + str(GP[2])
		id = str(rec.chrom)+":"+str(rec.pos) + ":" + rec.ref + ":" + rec.alts[0]
		gp_info[id] = GT + ";" + GP
		gp_info_GT[id] = GT

	gl_vcf_in = VariantFile(out + "/germline/" +  args.chr + ".gl.vcf.gz") 
	gl_vcf_dp4 = open(out + "/germline/" +  args.chr + ".gl.vcf.DP4","w")
	gl_vcf_filter_dp4 = open(out + "/germline/" +  args.chr + ".gl.vcf.filter.DP4","w")
	gl_vcf_filter_bed = open(out + "/germline/" +  args.chr + ".gl.vcf.filter.hc.bed","w")
	gl_vcf_filter_txt = open(out + "/germline/" +  args.chr + ".gl.vcf.filter.hc.pos","w")

	args.depth_filter_novelSNV = 10
	for rec in gl_vcf_in.fetch():
		info_var = rec.info['I16']
		id = str(rec.chrom)+":"+str(rec.pos) + ":" + rec.ref + ":" + rec.alts[0]
		gp_info_var = "NA"
		if (id in gp_info):
			gp_info_var = gp_info[id]
			base = rec.ref + "-" + rec.alts[0]
			if (gp_info_GT[id]=="0/0"): 
				allele_ratio = (info_var[2] + info_var[3])/(info_var[0] + info_var[1] + 1)
				if (base not in error_rate):
					error_rate[base] = 0
					error_rate_cnt[base] = 0
				error_rate[base]=allele_ratio + error_rate[base]
				error_rate_cnt[base] = error_rate_cnt[base] + 1

		a =  "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(rec.chrom, rec.pos, rec.ref, rec.alts[0], rec.info['DP'], info_var[0], info_var[1], info_var[2], info_var[3], gp_info_var)
		gl_vcf_dp4.write(a)

		# include all variants from germline 
		if (info_var[0] + info_var[1]>=4 and info_var[2] + info_var[3]>=4 or (id in gp_info)):
			tol_ref = info_var[0] + info_var[1]
			tol_alt = info_var[2] + info_var[3]		
			if(tol_ref > 0 and tol_alt/tol_ref>0.01):
				b = "{}\t{}\t{}\n".format(rec.chrom, rec.pos-1, rec.pos)
				gl_vcf_filter_bed.write(b)
				b = "{}\t{}\n".format(rec.chrom, rec.pos)
				gl_vcf_filter_txt.write(b)
				gl_vcf_filter_dp4.write(a)

def BamExtract(args):
	samtools = os.path.abspath(args.app_path) + "/samtools" 
	out = os.path.abspath(args.out)
	inbam = out + "/Bam/" + args.chr + ".filter.bam"
	outbam =  out + "/Bam/" + args.chr + ".filter.targeted.bam"
	# remembr to update bed file   ls -lrt

	cmd1 = samtools + " view " + " -b  -L " + out + "/germline/" +  args.chr + ".gl.vcf.filter.hc.bed " + inbam + " -o " + outbam
	with open(out+"/Script" + args.chr + "/BamExtract.sh","w") as f_out:
		f_out.write(cmd1 + "\n")
		f_out.write(samtools + " index " +  outbam + "\n")
	cmd="bash " + out+"/Script" + args.chr + "/BamExtract.sh"
	runCMD(cmd,args)



def robust_get_tag(read, tag_name):  
	try:  
		return read.get_tag(tag_name)
	except KeyError:
		return "NotFound"

def BamSplit(args):
	out = args.out
	bam_filter =  out + "/Bam/" + args.chr + ".filter.targeted.bam"
	assert os.path.isfile(bam_filter), "Bam filtering target file {} cannot be found! Please run germline mode firslty if you want to call somatic mutaitons".format(bam_filter)
	samtools = args.samtools 
	os.system("mkdir -p " + out + "/Bam/split_bam/")
	cell_clst = pd.read_csv(args.cell_cluster)   
	df = pd.DataFrame(cell_clst, columns= ['cell','cluster'])
	clst = df['cluster'].unique()
	clst_num = len(clst)
	bamlist = open(out + "/Bam/split_bam/" +  args.chr + ".file.lst","w")
	for i in range(clst_num):
		a = df['cell'][df['cluster']==clst[i]]
		#print(a)
		a = a.to_list()
		infile = pysam.AlignmentFile(bam_filter,"rb")
		# Note to change the read groups 

		tp =infile.header.to_dict()
		tp['RG'][0]['SM'] = "clst" + str(clst[i])
		tp['RG'][0]['ID'] = "clst" + str(clst[i])
		outfile =  pysam.AlignmentFile( out + "/Bam/split_bam/" + args.chr + "_" + str(clst[i]) + ".bam", "wb", header=tp)
		for s in infile:
			t  = robust_get_tag(s,"CB")
			if t in a:
				outfile.write(s)
		
		outfile.close()
		infile.close()
		cmd=samtools + " index " + out + "/Bam/split_bam/" + args.chr + "_" + str(clst[i]) + ".bam"
		runCMD(cmd,args)
		#os.system(samtools + " index " + out + "/Bam/split_bam/" + args.chr + "_" + str(clst[i]) + ".bam")
		bamlist.write(out + "/Bam/split_bam/" + args.chr + "_" + str(clst[i]) + ".bam" + "\n")

def runCMD(cmd):

	os.system(cmd)
	#process = subprocess.run(cmd, shell=True, stdout=open(args.logfile, 'w'), stderr=open(args.logfile,'w'))
	  

def validate_user_setting_somatic(args):

	assert os.path.isdir(args.out), "The germline output folder {} cannot be found! Please run germline module.".format(args.out)
	args.bam_filter = args.out + "/Bam/" + args.chr + ".filter.bam"
	assert os.path.isfile(args.bam_filter), "The filtered bam file {} cannot be found! Please run germline module".format(args.bam_filter)
	args.vcf_germline = args.out + "/germline/" + args.chr + ".germline.vcf"
	assert os.path.isfile(args.vcf_germline), "The germline vcf file {} cannot be found! Please run germline module".format(args.vcf_germline)
	args.vcf_depth = args.out + "/germline/" +  args.chr + ".gl.vcf.DP4"
	assert os.path.isfile(args.vcf_depth), "The sequencing depth file {} cannot be found! Please run germline module".format(args.vcf_depth)
	assert os.path.isfile(args.barcode), "The cell barcode file {} cannot be found!".format(args.barcode)


def vcf2mat(args):
	vcf_germline_phased = VariantFile(args.out + "/germline/" +  args.chr + "_phased.vcf.gz") 
	error  = open(args.out + "/germline/" +  args.chr + ".gl.vcf.filter.DP4", "r")
	vcf_in = VariantFile(args.out + "/germline/" +  args.chr + ".gl.filter.hc.cell.vcf.gz") 
	mat_out = gzip.open(args.out + "/germline/" +  args.chr + ".gl.filter.hc.cell.mat.gz","wt")
	allele_info = {}
	phase_info =  {}
	error_info ={}

	with open(args.out + "/germline/" +  args.chr + ".gl.vcf.filter.DP4",'r') as fp:
		for line in fp:
			data = line.split("\t")
			id=str(data[0])+":"+str(data[1])
			if re.search("0/0",data[9]):
				error_info[id] = "0|0"
				print(id)

	n = len(list((vcf_in.header.samples)))
	for rec in vcf_germline_phased.fetch():
		id = str(rec.chrom)+":"+str(rec.pos)
		allele =  rec.ref + ":" + rec.alts[0]
		phase = rec.samples.values()[0]['GT']
		phase = str(phase[0]) + "|" + str(phase[1])
		allele_info[id] = allele
		phase_info[id] = phase

	for rec in vcf_in.fetch():
		#print(rec.samples.values())

		id = str(rec.chrom)+":"+str(rec.pos)
		allele =  rec.ref + ":" + rec.alts[0]
		DP = rec.info["DP"]
		a = [None]*(4+n)
		a[0] = id
		a[1] = allele
		a[2] = str(DP)
		a[3] = ".|."
		matched_allele = 1 
		i = 3

		# note the 1|1 genotypes are not used 
		if (id in phase_info):
			a[3] = phase_info[id]
			if not allele==allele_info[id]:
				matched_allele == 0
			if matched_allele and phase_info[id]=="0|1": 
				for value in rec.samples.values():
					ref = value['DP4'][0] +  value['DP4'][1]
					alt = value['DP4'][2] +  value['DP4'][3]
					i = i + 1
					a[i] = str(ref)+"|"+str(alt)
				mat_out.write('\t'.join(a))
				mat_out.write('\n')
			if matched_allele and phase_info[id]=="1|0": 
				for value in rec.samples.values():
					ref = value['DP4'][0] +  value['DP4'][1]
					alt = value['DP4'][2] +  value['DP4'][3]
					i = i + 1
					a[i] = str(alt)+"|"+str(ref)
				mat_out.write('\t'.join(a))
				mat_out.write('\n')
		else:
			if (id in error_info):
				a[3] = error_info[id]
			for value in rec.samples.values():
				ref = value['DP4'][0] +  value['DP4'][1]
				alt = value['DP4'][2] +  value['DP4'][3]
				i = i + 1
				a[i] = str(ref)+"/"+str(alt)
			mat_out.write('\t'.join(a))
			mat_out.write('\n')

	mat_out.close()
	vcf_in.close()




def bam2mat(args): 

	samtools = os.path.abspath(args.app_path) + "/samtools"
	beagle =  os.path.abspath(args.app_path) + "/beagle.27Jul16.86a.jar" 
	bam_filter =  args.out + "/Bam/" + args.chr + ".filter.targeted.bam"
	snv_pos =  args.out + "/germline/" + args.chr + ".gl.vcf.filter.hc.bed"
	bam_lst = args.out + "/Bam/cell_bam.lst"
	vcf_out = args.out + "/germline/" + args.chr + ".gl.filter.hc.cell.vcf.gz"
	out = args.out

	cmd = samtools + " mpileup  -u -q 20 -Q 20  -t DP4   -d 10000000  -l " + snv_pos + " -b " + bam_lst + " -f  /rsrch3/scratch/bcb/jdou1/scAncestry/ref/fasta/genome.fa | " +  args.bcftools + " view | bgzip -c > " + vcf_out
	args.map = "/rsrch3/scratch/bcb/jdou1/scAncestry/ref/1KG3/plink." + args.chr + ".phase.addchr.GRCh38.map"
	args.imputation_panel  = "/rsrch3/scratch/bcb/jdou1/scAncestry/ref/1KG3/CCDG_14151_B01_GRM_WGS_2020-08-05_" + args.chr + ".filtered.shapeit2-duohmm-phased.vcf.gz"
	cmd3 = args.java + " -Xmx20g -jar " + beagle +  " gt=" +  out + "/germline/" +  args.chr + ".gt.vcf.gz  map="  + args.map +  " ref=" +  args.imputation_panel  + "  chrom=" + args.chr  + " out="   +  args.out + "/germline/" + args.chr + "_phased " + "impute=false  modelscale=2  nthreads=48  gprobs=true  niterations=0"
		
	with open(args.out+"/Script" + args.chr + "/Bam2mat.sh","w") as f_out:
		#f_out.write(cmd3 + "\n")
		f_out.write(cmd + "\n")
	cmd="bash " + args.out+"/Script" + args.chr + "/Bam2mat.sh"
	
	runCMD(cmd,args)
	vcf2mat(args)

    
#def somatic_haplotype(args):


