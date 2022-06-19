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




def print_parameters_given(args):
	logger.info("Parameters in effect:")
	for arg in vars(args):
		if arg=="func": continue
		logger.info("--{} = [{}]".format(arg, vars(args)[arg]))


def get_seq_type_from_user_cfg(fn_user_cfg):
	with open(fn_user_cfg) as fh:
		user_cfg = dict(yaml.safe_load(fh))
	seq_type=user_cfg["seqType"]
	if (seq_type!="WES") and (seq_type!="WGS"):
		logger.error("seqType can only be WES or WGS! exiting!")
		exit(1)

	return seq_type


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


def validate_user_setting(args):
	assert os.path.isfile(args.bamFile), "The bam list file {} cannot be found!".format(args.bam_list)
	assert os.path.isfile(args.reference), "The genome reference fasta file {} cannot be found!".format(args.reference)
	assert os.path.isfile(args.imputation_panel), "Filtered genotype file of 1KG3 ref panel {} cannot be found!".format(args.imputation_panel)
	

def check_dependencies(args):
	programs_to_check = ("vcftools", "bcftools", "beagle.08Feb22.fa4.jar", "beagle.27Jul16.86a.jar","samtools","picard.jar", "java")

	for prog in programs_to_check:
		out = os.popen("command -v {}".format(args.app_path + "/" + prog)).read()
		assert out != "", "Program {} cannot be found!".format(prog)

#	python_pkgs_to_check = ("drmaa",)

#	for pkg in python_pkgs_to_check:
#		out_pipe = os.popen('python -c "import {}"'.format(pkg))

#		assert out_pipe.close() is None, "Python module {} has not been installed!".format(pkg)


def BamFilter(args):
		
	infile = pysam.AlignmentFile(args.bamFile,"rb")
	outfile =  pysam.AlignmentFile( args.out + "/Bam/" + args.chr + ".filter.bam", "wb", template=infile)
	for s in infile.fetch(args.chr):
		if (s.query_length - s.get_tag("AS")) < args.max_softClipped  and s.get_tag("nM") < args.max_mismatch:
			outfile.write(s)
	infile.close()
	outfile.close()
	os.system(args.samtools + " index " +  args.out + "/Bam/" + args.chr + ".filter.bam")


def getDPinfo(args):
	out = args.out
	gp_vcf_in = VariantFile(out + "/SCvarCall/" +  args.chr + ".gp.vcf.gz") 
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

	gl_vcf_in = VariantFile(out + "/SCvarCall/" +  args.chr + ".gl.vcf.gz") 
	gl_vcf_dp4 = open(out + "/SCvarCall/" +  args.chr + ".gl.vcf.DP4","w")
	gl_vcf_filter_dp4 = open(out + "/SCvarCall/" +  args.chr + ".gl.vcf.filter.DP4","w")
	gl_vcf_filter_bed = open(out + "/SCvarCall/" +  args.chr + ".gl.vcf.filter.hc.bed","w")
	gl_vcf_filter_txt = open(out + "/SCvarCall/" +  args.chr + ".gl.vcf.filter.hc.pos","w")

	for rec in gl_vcf_in.fetch():
		info_var = rec.info['I16']
		id = str(rec.chrom)+":"+str(rec.pos) + ":" + rec.ref + ":" + rec.alts[0]
		gp_info_var = "NA"
		if (id in gp_info):
			gp_info_var = gp_info[id]
			base = rec.ref + "-" + rec.alts[0]
			if (gp_info_GT[id]=="0/0"): 
				allele_ratio = (info_var[2] + info_var[3])/(info_var[0] + info_var[1] + 1)
				#print(id + " " + gp_info[id] + base)
				if (base not in error_rate):
					error_rate[base] = 0
					error_rate_cnt[base] = 0
				error_rate[base]=allele_ratio + error_rate[base]
				error_rate_cnt[base] = error_rate_cnt[base] + 1

		a =  "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(rec.chrom, rec.pos, rec.ref, rec.alts[0], rec.info['DP'], info_var[0], info_var[1], info_var[2], info_var[3], gp_info_var)
		gl_vcf_dp4.write(a)
		if (info_var[0] + info_var[1]>=4 and info_var[2] + info_var[3]>=4):
			gl_vcf_filter_dp4.write(a)
			tol_ref = info_var[0] + info_var[1]
			tol_alt = info_var[2] + info_var[3]
			
			if(tol_ref > args.depth_filter and tol_alt/tol_ref>args.alt_ratio):
				b = "{}\t{}\t{}\n".format(rec.chrom, rec.pos-200, rec.pos+200)
				gl_vcf_filter_bed.write(b)
				b = "{}\t{}\n".format(rec.chrom, rec.pos)
				gl_vcf_filter_txt.write(b)

def BamExtract(args):
	samtools = os.path.abspath(args.app_path) + "/samtools" 
	out = os.path.abspath(args.out)
	inbam = out + "/Bam/" + args.chr + ".filter.bam"
	outbam =  out + "/Bam/" + args.chr + ".filter.targeted.bam"
	# remembr to update bed file   ls -lrt

	cmd1 = samtools + " view " + " -b  -L " + out + "/SCvarCall/" +  args.chr + ".gl.vcf.filter.hc.bed " + inbam + " -o " + outbam
	with open(out+"/Script/ExtractBam.sh","w") as f_out:
		f_out.write(cmd1 + "\n")
		f_out.write(samtools + " index " +  outbam + "\n")
	cmd="bash " + out+"/Script/ExtractBam.sh"
	runCMD(cmd,args)

	#for p in error_rate:
	#	print(p + "---" + str(error_rate[p]) + "---" + str(error_rate_cnt[p]))


def robust_get_tag(read, tag_name):  
	try:  
		return read.get_tag(tag_name)
	except KeyError:
		return "NotFound"

def BamSplit(args):
	out = args.out
	bam_filter =  out + "/Bam/" + args.chr + ".filter.targeted.bam"
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

def RunMonova(args):

	# run Monova for cluster-level mutation calling 
	samtools = args.samtools 
	snp_pos = args.out + "/SCvarCall/" +  args.chr + ".gl.vcf.filter.hc.pos"
	bamlst = args.out +  "/Bam/split_bam/" +  args.chr + ".file.lst"
	ref = args.reference
	monovar =  os.path.abspath(args.app_path) + "/../src/monovar.py"
	outvcf = args.out + "/SCvarCall/" +  args.chr + ".monovar.vcf"
	cmd = samtools + " mpileup   -BQ0 -d10000 -q 40   -l " + snp_pos +  " -f " + ref + " -b " + bamlst + " | python " +   monovar +   " -f " + ref + " -b " + bamlst +" -p 0.002 -a 0.2 -t 0.05 -m 2  -o " + outvcf
	with open(args.out+"/Script/runMonova.sh","w") as f_out:
			f_out.write(cmd + "\n")
	cmd = "bash " + args.out+"/Script/runMonova.sh"
	runCMD(cmd,args)

def runCMD(cmd, args):
	#print(cmd)
	os.system(cmd + " > " + args.logfile)
	#process = subprocess.run(cmd, shell=True, stdout=open(args.logfile, 'w'), stderr=open(args.logfile,'w'))
	

def SCvarCall(args):

	out = args.out
	logger.info("Preparing varint calling pipeline...")
	print_parameters_given(args)

	logger.info("Checking existence of essenstial resource files...")
	validate_user_setting(args)

	logger.info("Checking dependencies...")
	check_dependencies(args)

	
	os.system("mkdir -p " + out )
	os.system("mkdir -p " + out +  "/Bam")
	os.system("mkdir -p " + out +  "/SCvarCall")
	os.system("mkdir -p " + out +  "/Script")
	
	samtools = args.samtools 
	bcftools = args.bcftools 
	#java = os.path.abspath(args.app_path) + "/java"
	java =  args.java
	beagle = args.beagle
	
	logger.info("Filtering bam files...")
	BamFilter(args)

	cmd1 = samtools + " mpileup " + args.bamFile + " -f "  + args.reference  + " -r " +  args.chr + " -q 20 -Q 20 -t DP -d 10000000 -v "
	cmd1 = cmd1 + " | " + bcftools + " view " + " | "  + bcftools  + " norm -m-both -f " + args.reference 
	cmd1 = cmd1 + " | grep -v X | grep -v INDEL | bgzip -c > " + out + "/SCvarCall/" +  args.chr + ".gl.vcf.gz" 

	cmd2 = bcftools + " view " +  out + "/SCvarCall/" +  args.chr + ".gl.vcf.gz" + " -i 'FORMAT/DP>10' | " + bcftools + " call -cv  | bgzip -c > " +  out + "/SCvarCall/"  +  args.chr + ".gt.vcf.gz"
	cmd3 = java + " -Xmx20g -jar " + beagle +  " gl=" +  out + "/SCvarCall/" +  args.chr + ".gl.vcf.gz"  +  " ref=" +  args.imputation_panel  + "  chrom=" + args.chr  + " out="   +  out + "/SCvarCall/" + args.chr + ".gp " + "impute=false  modelscale=2  nthreads=48  gprobs=true  niterations=0"
	cmd4 = "zless -S " +  out + "/SCvarCall/" + args.chr + ".gp.vcf.gz | grep -v  0/0  > " +  out + "/SCvarCall/" + args.chr + ".germline.vcf"

	with open(out+"/Script/runBeagle.sh","w") as f_out:
			f_out.write(cmd1 + "\n")
			f_out.write(cmd2 + "\n")
			f_out.write(cmd3 + "\n")
			f_out.write(cmd4 + "\n")

	logger.info("Performing Variant Calling...")
	cmd = "bash " + out+"/Script/runBeagle.sh"
	runCMD(cmd,args)
	logger.info("Generating variant statistical information ...")
	getDPinfo(args)
	logger.info("Extracting reads in potential somatic regions ...")
	BamExtract(args)
	logger.info("Splitting reads baed on cell cluster ID ...")
	BamSplit(args)
	logger.info("Running monova ...")
	RunMonova(args)



def main():
	parser = argparse.ArgumentParser(
		description="""scPopGene: single-cell population genetics analysis)
		""",
		epilog=
		"""Typical workflow: SCvarCall => SCancestry => SCqtl\nCitation: TBD
		""",
		formatter_class=argparse.RawTextHelpFormatter)

	subparsers = parser.add_subparsers(title='Available subcommands', dest="subcommand")
	
	# every subcommand needs user config file
	common_parser = argparse.ArgumentParser(add_help=False)
	parser_varCall = subparsers.add_parser('SCvarCall', parents=[common_parser],
		help='Variant discovery, genotype calling from single cell data (scRNA-seq, snRNA-seq or snATAC-seq)',
		formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser_varCall.add_argument('-b', '--bamFile', required=True,
								help="The bam file for the study sample, the bam file should be sorted")
	parser_varCall.add_argument('-c', '--chr', required= True, 
								help="The chromosome used for variant calling")
	parser_varCall.add_argument('-o', '--out', required= False,
								help="The output director")
	parser_varCall.add_argument('-r', '--reference', required= True, 
								help="The human genome reference used for alignment")
	parser_varCall.add_argument('-p', '--imputation-panel', required= True, 
								help="The population-level variant panel for variant refinement such as 1000 Genome 3")
	parser_varCall.add_argument('-d', '--depth_filter', required=False, type=int, default=50,
								help="The sequencing depth filter for variants not overlapped with public database")
	parser_varCall.add_argument('-t', '--alt_ratio', required=False, type=float, default=0.1,
								help="The minina allele frequency for variants as potential somatic mutation")
	parser_varCall.add_argument('-m', '--max-mismatch', required=False, type=int, default=3,
								help="The maximal mismatch allowed in one reads for variant calling")
	parser_varCall.add_argument('-s', '--max-softClipped', required=False, type=int, default=1,
								help="The maximal soft-clipped allowed in one reads for variant calling")
	parser_varCall.add_argument('-a', '--app-path', required=True,
								help="The app library paths used in the tool")
	parser_varCall.add_argument('-i', '--cell_cluster', required=True,
								help="The cell cluster csv file used for somatic variant calling")
	parser_varCall.set_defaults(func=SCvarCall)

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
	args.java =  "java"
	args.beagle = os.path.abspath(args.app_path) + "/beagle.27Jul16.86a.jar"
	args.bamFile = os.path.abspath(args.bamFile)
	args.func(args)

	logger.info("Success! See instructions above.")

if __name__ == "__main__":
	main()