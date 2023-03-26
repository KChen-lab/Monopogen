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
	assert os.path.isfile(args.imputation_panel), "Filtered genotype file of 1KG3 ref panel {} cannot be found!".format(args.imputation_panel)
	if args.mode=="single":
		assert os.path.isfile(args.bamFile), "The bam file {} cannot be found!".format(args.bamFile)
		assert os.path.isfile(args.bamFile + ".bai"), "The bam.bai file {} cannot be found!".format(args.bamFile)
	if args.mode=="multiple":
		assert os.path.isfile(args.bamFile), "The bam file {} cannot be found!".format(args.bamFile)
		# check whether each bam file available	
		with open(args.bamFile) as fh:
			assert os.path.isfile(fh), "The bam file {} cannot be found!".format(fh)
			assert os.path.isfile(fh + ".bai"), "The bam.bai file {} cannot be found!".format(fh)


def check_dependencies(args):
	programs_to_check = ("vcftools", "bgzip",  "bcftools", "beagle.08Feb22.fa4.jar", "beagle.27Jul16.86a.jar","samtools","picard.jar", "java")

	for prog in programs_to_check:
		out = os.popen("command -v {}".format(args.app_path + "/" + prog)).read()
		assert out != "", "Program {} cannot be found!".format(prog)

#	python_pkgs_to_check = ("drmaa",)

#	for pkg in python_pkgs_to_check:
#		out_pipe = os.popen('python -c "import {}"'.format(pkg))

#		assert out_pipe.close() is None, "Python module {} has not been installed!".format(pkg)



def addChr(in_bam):
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
	os.system(args.samtools + " index " +  out_bam)
	os.system(" mv " + out_bam + " " + in_bam)
	os.system(" mv " + out_bam + ".bai  " + in_bam + ".bai")



def BamFilter(args):
	infile = pysam.AlignmentFile(args.bamFile,"rb")
	contig_names = infile.references
	cnt=0 
	search_chr = args.chr
	for contig in contig_names:
		if contig.startswith("chr"):
			cnt=cnt+1
	if cnt==0:
			logger.info("The contig {} does not contain the prefix 'chr' and we will add 'chr' on it ".format(args.bamFile))
			search_chr = args.chr[3:]
			#searchBam = args.bamFile+".bam"
			#addChr(args.bamFile, searchBam)
			#infile.close()
			#infile = pysam.AlignmentFile(searchBam,"rb")

	tp =infile.header.to_dict()
			
	#print(tp)

	if not "RG" in tp:
		sampleID = os.path.splitext(os.path.basename(args.bamFile))[0]
		tp1 = [{'SM':sampleID,'ID':sampleID, 'LB':"0.1", 'PL':"ILLUMINA", 'PU':sampleID}]
		tp.update({'RG': tp1})
		#print(tp)
		#tp['RG'][0]['SM'] = 
		#tp['RG'][0]['ID'] = os.path.splitext(os.path.basename(args.bamFile))[0]
	#print(tp)   
	outfile =  pysam.AlignmentFile( args.out + "/Bam/" + args.chr + ".filter.bam", "wb", header=tp)


	for s in infile.fetch(search_chr):  
		#print(str(s.query_length)  + ":" + str(s.get_tag("AS")) + ":" + str(s.get_tag("NM")))
		if s.has_tag("NM"):
			val= s.get_tag("NM")
		if s.has_tag("nM"):
			val= s.get_tag("nM")                  
		if val < args.max_mismatch:
			outfile.write(s)
	infile.close()
	outfile.close()

	os.system(args.samtools + " index " +  args.out + "/Bam/" + args.chr + ".filter.bam")
	if cnt ==0:
		addChr(args.out + "/Bam/" + args.chr + ".filter.bam")
	args.bam_filter = args.out + "/Bam/" + args.chr + ".filter.bam"

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

def runCMD(cmd, args):
	#print(cmd)
	os.system(cmd + " > " + args.logfile)
	#process = subprocess.run(cmd, shell=True, stdout=open(args.logfile, 'w'), stderr=open(args.logfile,'w'))
	  

def germline(args):
	logger.info("Performing germline variant calling...")
	print_parameters_given(args)

	logger.info("Checking existence of essenstial resource files...")
	validate_user_setting_germline(args)

	logger.info("Checking dependencies...")
	check_dependencies(args)
	out = args.out
	os.system("mkdir -p " + out )
	os.system("mkdir -p " + out +  "/Bam")
	os.system("mkdir -p " + out +  "/germline")
	os.system("mkdir -p " + out +  "/Script" + args.chr)
	
	samtools = args.samtools 
	bcftools = args.bcftools 
	#java = os.path.abspath(args.app_path) + "/java"
	java =  args.java
	beagle = args.beagle
	if args.mode == "single":
		if args.step == "bamQC" or args.step == "all":
			logger.info("Filtering bam files...")
			BamFilter(args)   
		
		args.bam_filter = args.out + "/Bam/" + args.chr + ".filter.bam"

		cmd1 = samtools + " mpileup " + args.bam_filter + " -f "  + args.reference  + " -r " +  args.chr + " -q 20 -Q 20 -t DP -d 10000000 -v "
		cmd1 = cmd1 + " | " + bcftools + " view " + " | "  + bcftools  + " norm -m-both -f " + args.reference 
		cmd1 = cmd1 + " | grep -v \"<X>\" | grep -v INDEL |" + args.bgzip +   " -c > " + out + "/germline/" +  args.chr + ".gl.vcf.gz" 
		cmd2 = bcftools + " view " +  out + "/germline/" +  args.chr + ".gl.vcf.gz" + " -i 'FORMAT/DP>1' | " + bcftools + " call -cv  | " + args.bgzip +    "  -c > " +  out + "/SCvarCall/"  +  args.chr + ".gt.vcf.gz"
		cmd3 = java + " -Xmx20g -jar " + beagle +  " gl=" +  out + "/germline/" +  args.chr + ".gl.vcf.gz"  +  " ref=" +  args.imputation_panel  + "  chrom=" + args.chr  + " out="   +  out + "/germline/" + args.chr + ".gp " + "impute=false  modelscale=2  nthreads=48  gprobs=true  niterations=0"
		
		cmd4 = "zless -S " +  out + "/germline/" + args.chr + ".gp.vcf.gz | grep -v  0/0  > " +  out + "/germline/" + args.chr + ".germline.vcf"
		cmd5 = java + " -Xmx20g -jar " + beagle +  " gt=" +  out + "/germline/" +  args.chr + ".germline.vcf"  +  " ref=" +  args.imputation_panel   +  "  chrom=" + args.chr  + " out="   +  out + "/germline/" + args.chr + ".phased " + "impute=false  modelscale=2  nthreads=48  gprobs=true  niterations=0"
		

		if args.step == "varScan" or args.step == "all":
			with open(out+"/Script" + args.chr + "/runVarScan.sh","w") as f_out:
				#if args.step=="germline" or args.step=="all":
				f_out.write(cmd1 + "\n")
				f_out.write(cmd2 + "\n")
			cmd = "bash " + out+"/Script" + args.chr +  "/runVarScan.sh"
			runCMD(cmd,args)

		if args.step == "varImpute" or args.step == "all":
			with open(out+"/Script" + args.chr + "/runVarImpute.sh","w") as f_out:
				f_out.write(cmd3 + "\n")
				f_out.write(cmd4 + "\n")
			cmd = "bash " + out+"/Script" + args.chr +  "/runVarImpute.sh"
			print(cmd)
			runCMD(cmd,args)
			getDPinfo(args)
		if args.step == "varPhasing" or args.step == "all":
			with open(out+"/Script" + args.chr + "/runVarPhasing.sh","w") as f_out:
				f_out.write(cmd5 + "\n")
			cmd = "bash " + out+"/Script" + args.chr +  "/runVarPhasing.sh"
			print(cmd)
			runCMD(cmd,args)
			getDPinfo(args)

      	
	if args.mode == "multi":
		#if args.step == "bamQC" or args.step == "all":
		#	logger.info("Filtering bam files...")
			#BamFilter(args.bamFile)
		if args.step == "varScan" or args.step == "all" or args.stepp == "varImpute":
			cmd1 = samtools + " mpileup -b " + args.bamFile + " -f "  + args.reference  + " -r " +  args.chr + " -q 20 -Q 20 -t DP -d 10000000 -v "
			cmd1 = cmd1 + " | " + bcftools + " view " + " | "  + bcftools  + " norm -m-both -f " + args.reference 
			cmd1 = cmd1 + " | grep -v \"<X>\"  | grep -v INDEL |" + args.bgzip +   " -c > " + out + "/germline/" +  args.chr + ".gl.vcf.gz" 
			cmd2 = bcftools + " view " +  out + "/germline/" +  args.chr + ".gl.vcf.gz" + " -i 'FORMAT/DP>10' | " + bcftools + " call -cv  | " + args.bgzip +    "  -c > " +  out + "/SCvarCall/"  +  args.chr + ".gt.vcf.gz"
			cmd3 = java + " -Xmx20g -jar " + beagle +  " gl=" +  out + "/germline/" +  args.chr + ".gl.vcf.gz"  +  " ref=" +  args.imputation_panel  + "  chrom=" + args.chr  + " out="   +  out + "/germline/" + args.chr + ".germline " + "impute=false  modelscale=2  nthreads=48  gprobs=true  niterations=0"
			#cmd4 = "zless -S " +  out + "/SCvarCall/" + args.chr + ".gp.vcf.gz | grep -v  0/0  > " +  out + "/SCvarCall/" + args.chr + ".germline.vcf"

			with open(out+"/Script" + args.chr + "runVarScan.sh","w") as f_out:
					f_out.write(cmd1 + "\n")
					f_out.write(cmd3 + "\n")
			cmd = "bash " + out+"/Script" + args.chr +  "/runVarScan.sh"
			runCMD(cmd,args)


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



def somatic(args):
	
	validate_user_setting_somatic(args)
	getDPinfo(args)
	BamExtract(args)
	# assume that we have split bam file into cell level
	bam2mat(args)






def main():
	parser = argparse.ArgumentParser(
		description="""Monopogen: SNV calling from single cell sequencing
		""",
		epilog=
		"""Typical modules: germline, somatic
		""",
		formatter_class=argparse.RawTextHelpFormatter)

	subparsers = parser.add_subparsers(title='Available subcommands', dest="subcommand")
	
	# every subcommand needs user config file
	common_parser = argparse.ArgumentParser(add_help=False)
	parser_germline = subparsers.add_parser('germline', parents=[common_parser],
		help='Germline variant discovery, genotype calling from single cell data',
		formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser_germline.add_argument('-b', '--bamFile', required=True,
								help="The bam file for the study sample, the bam file should be sorted") 
	parser_germline.add_argument('-y', '--mode', required=True,
								choices=['single','multi'],
								help="Single sample or multiple samples. Only available for germline variant calling mode. This step can increase variant detection.")
	parser_germline.add_argument('-c', '--chr', required= True, 
								help="The chromosome used for variant calling")
	parser_germline.add_argument('-t', '--step', required= True, default="all",
								choices=['bamQC', 'varScan', 'varImpute' , 'varPhasing', 'all'],
								help="Run germline variant calling step by step")
	parser_germline.add_argument('-o', '--out', required= False,
								help="The output director")
	parser_germline.add_argument('-r', '--reference', required= True, 
								help="The human genome reference used for alignment")
	parser_germline.add_argument('-p', '--imputation-panel', required= True, 
								help="The population-level variant panel for variant imputation refinement, such as 1000 Genome 3")
	parser_germline.add_argument('-d', '--depth_filter_novelSNV', required=False, type=int, default=24,
								help="The minimal read depth supported to call novel SNVs not listed in reference panel")
	parser_germline.add_argument('-m', '--max-mismatch', required=False, type=int, default=3,
								help="The maximal mismatch allowed in one reads for variant calling")
	parser_germline.add_argument('-s', '--max-softClipped', required=False, type=int, default=1,
								help="The maximal soft-clipped allowed in one reads for variant calling")
	parser_germline.add_argument('-a', '--app-path', required=True,
								help="The app library paths used in the tool")
	parser_germline.set_defaults(func=germline)



	parser_somatic = subparsers.add_parser('somatic', parents=[common_parser],
		help='Somatic variant calling from single cell data',
		formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser_somatic.add_argument('-i', '--input-folder', required=True,
								help="The output folder from previous germline module")
	parser_somatic.add_argument('-c', '--chr', required= True, 
								help="The chromosome used for variant calling")
	parser_somatic.add_argument('-l', '--barcode', required= True, 
								help="The csv file including cell barcode information")
	parser_somatic.add_argument('-a', '--app-path', required=True,
								help="The app library paths used in the tool")
	parser_somatic.set_defaults(func=somatic)

	args = parser.parse_args()
	if args.subcommand is None:
		# if no command is specified, print help and exit
		print("Please specify one subcommand! Exiting!")
		print("-"*80)
		parser.print_help()
		exit(1)

	# check if user configure exists
	#assert os.path.isfile(args.userCfg), ("User config file {} cannot be found".format(args.userCfg))

	# execute subcommand-specific function
	print(args.subcommand)

	if args.subcommand == "somatic":
		args.out = args.input_folder

	args.logfile = args.out + "_" + args.chr + ".log"
	if os.path.exists(args.logfile):
		os.remove(args.logfile)
	handler1 = logging.FileHandler(args.logfile)
	handler1.setFormatter(logging.Formatter(
		'[{asctime}] {levelname:8s} {filename} {message}', style='{'))
	logger.addHandler(handler1)

	args.out = os.path.abspath(args.out)
	args.samtools  = os.path.abspath(args.app_path) + "/samtools" 
	args.bcftools = os.path.abspath(args.app_path) + "/bcftools"
	args.bgzip = os.path.abspath(args.app_path) + "/bgzip"
	args.java =  "java"
	args.beagle = os.path.abspath(args.app_path) + "/beagle.27Jul16.86a.jar"
	if args.subcommand == "germline":
		args.bamFile = os.path.abspath(args.bamFile)
	args.func(args)

	logger.info("Success! See instructions above.")








if __name__ == "__main__":
	main()
