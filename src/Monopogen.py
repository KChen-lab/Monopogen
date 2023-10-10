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
from germline import *
from somatic import *
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




def germline(args):
	logger.info("Performing germline variant calling...")
	print_parameters_given(args)

	logger.info("Checking existence of essenstial resource files...")
	validate_user_setting_germline(args)

	logger.info("Checking dependencies...")
	check_dependencies(args)
	out = args.out
	os.system("mkdir -p " + out )
	
	os.system("mkdir -p " + out +  "/germline")
	os.system("mkdir -p " + out +  "/Script")


	# check whether region files were set correctly 
	joblst = []
	with open(args.region) as f_in:
		for line in f_in:
			record = line.strip().split(",")
			if(len(record)==1):
				jobid = record[0]
			if(len(record)==3):
				jobid = record[0] + ":" + record[1] + "-" + record[2]
			bam_filter = args.out + "/Bam/" +  record[0] +  ".filter.bam.lst"
			imputation_vcf = args.imputation_panel + "CCDG_14151_B01_GRM_WGS_2020-08-05_" + record[0] + ".filtered.shapeit2-duohmm-phased.vcf.gz"
			cmd1 = samtools + " mpileup -b " + bam_filter + " -f "  + args.reference  + " -r " +  jobid + " -q 20 -Q 20 -t DP -d 10000000 -v "
			cmd1 = cmd1 + " | " + bcftools + " view " + " | "  + bcftools  + " norm -m-both -f " + args.reference 
			cmd1 = cmd1 + " | grep -v \"<X>\" | grep -v INDEL |" + bgzip +   " -c > " + args.out + "/germline/" +  jobid + ".gl.vcf.gz" 
			#cmd2 = bcftools + " view " +  out + "/germline/" +  jobid + ".gl.vcf.gz" + " -i 'FORMAT/DP>1' | " + bcftools + " call -cv  | " + bgzip +    "  -c > " +  args.out + "/SCvarCall/"  +  jobid + ".gt.vcf.gz"
			cmd3 = java + " -Xmx20g -jar " + beagle +  " gl=" +  out + "/germline/" +  jobid + ".gl.vcf.gz"  +  " ref=" +  imputation_vcf   + "  chrom=" + record[0] + " out="   +  out + "/germline/" + jobid + ".gp " + "impute=false  modelscale=2  nthreads=24  gprobs=true  niterations=0"
			
			cmd4 = "zless -S " +  out + "/germline/" + jobid + ".gp.vcf.gz | grep -v  0/0  > " +  out + "/germline/" + jobid + ".germline.vcf"
			cmd5 = java + " -Xmx20g -jar " + beagle +  " gt=" +  out + "/germline/" +  jobid + ".germline.vcf"  +  " ref=" +  imputation_vcf    +  "  chrom=" + record[0]  + " out="   +  out + "/germline/" + jobid+ ".phased " + "impute=false  modelscale=2  nthreads=24  gprobs=true  niterations=0"
			f_out = open(out + "/Script/runGermline_" +  jobid +  ".sh","w")
			if args.step == "varScan" or args.step == "all":
				f_out.write(cmd1 + "\n")
			#NSNV = withSNVs(out + "/germline/" +  jobid + ".gl.vcf.gz")
				#f_out.write(cmd2 + "\n")
			if args.step == "varImpute" or args.step == "all":
				#if NSNV>100:
					f_out.write(cmd3 + "\n")
					f_out.write(cmd4 + "\n")
			if args.step == "varPhasing" or args.step == "all":
				#if NSNV>100:
					f_out.write(cmd5 + "\n")
			
			joblst.append("bash " + out + "/Script/runGermline_" +  jobid +  ".sh")
	f_out.close()

	if not args.norun == "TRUE":
		with Pool(processes=args.nthreads) as pool:
			print(joblst)
			result = pool.map(runCMD, joblst)
	#error_check(all = region_lst, output = result, step = "germline module")



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


def somatic(args):
	
	validate_user_setting_somatic(args)
	os.system("mkdir -p " + args.out +  "/somatic")

	chr_lst = []
	region_lst = []
	with open(args.region) as f_in:
		for line in f_in:
			record = line.strip().split(",")
			if(len(record)==1):
			 	region = record[0]
			if(len(record)==3):
				region = record[0] + ":" + record[1] + "-" + record[2]
			chr_lst.append(record[0])
			region_lst.append(region)
	chr_lst = list(set(chr_lst))

	

	if args.step=="featureInfo" or args.step=="all":
		logger.info("Get feature information from sequencing data...")
		
		joblst = []
		for id in region_lst:
			joblst.append(id+">"+args.out)
		with Pool(processes=args.nthreads) as pool:
			result = pool.map(featureInfo, joblst)
		
	if args.step=="cellScan" or args.step=="all":
		
		logger.info("Get single cell level information from sequencing data...")
		
		chr_lst = sort_chr(chr_lst)

		run = 1
		if run:
			joblst = []
			for id in chr_lst:
				joblst.append(id+">"+args.out+">"+args.app_path)
			with Pool(processes=args.nthreads) as pool:
				result = pool.map(bamExtract, joblst)

			####### merge bams from different chromosomes 
			bamlst = []
			print(chr_lst)
			for chr in chr_lst:
				bam_filter =  out + "/Bam/" + chr + ".filter.targeted.bam"
				bamlst.append(bam_filter)
			print(bamlst)
			output_bam = out + "/Bam/merge.filter.targeted.bam"
			pysam.merge("-f","-o",output_bam,*bamlst)
			os.system(samtools + " index " + output_bam)

		if run:
			os.system("mkdir -p " + out + "/Bam/split_bam/")
			cell_clst = pd.read_csv(args.barcode)   
			df = pd.DataFrame(cell_clst, columns= ['cell','id'])
			df = df.sort_values(by=['id'])
			args.keep = float(args.keep)
			if args.keep < 1:
				dis = np.cumsum(df['id'])/np.sum(df['id'])
				N = sum(dis>(1-args.keep))
				df = df.iloc[-(N):]
			cell_lst = df['cell'].unique()
			joblst = []

			for cell in cell_lst:
					para = "merge" + ":" + cell + ":" + args.out + ":" + args.app_path
					joblst.append(para)

			with Pool(processes=args.nthreads) as pool:
				result = pool.map(bamSplit, joblst)

			if sum(result)==0:
				logger.error("No reads detected for cells in " + args.barcode + ". Please check 1) the input cell barcode file is matched with bam file; 2) the cell barcode name has the same format shown in bam file. For example XX-1!")
				logger.error("Failed! See instructions above.")
				exit(1)

			# generate the bam file list 
			cell_bam = open(out + "/Bam/split_bam/cell.bam.lst","w")
			for cell in cell_lst:
				cell_bam.write(out + "/Bam/split_bam/" + cell + ".bam\n")
			cell_bam.close()


		region_lst = []
		if args.winSize=="10MB":
			region_file = args.app_path + "/../resource/GRCh38.region.10MB.lst"
		if args.winSize=="50MB":
			region_file = args.app_path + "/../resource/GRCh38.region.50MB.lst"
		with open(region_file) as f_in:
			for line in f_in:
				record = line.strip().split(",")
				if(len(record)==1):
				 	region = record[0]
				if(len(record)==3):
					region = record[0] + ":" + record[1] + "-" + record[2]
				if record[0] in chr_lst:
					region_lst.append(region)

		joblst = []
		for id in region_lst:
			record = id.strip().split(":")
			chr = record[0]
			joblst.append(id+">"+chr+">"+args.out+">"+args.app_path+">"+args.reference)


		with Pool(processes=args.nthreads) as pool:
			result = pool.map(jointCall, joblst)
		error_check(all = region_lst, output = result, step = "cellScan:joint calling")


		##### merge vcfs from multiple regions 
	
		for id in chr_lst:
			tp = ""
			for reg in region_lst:
				if re.search(id+":", reg):
					tp = tp + " " + args.out+"/somatic/" +  reg + ".cell.gl.vcf.gz"
			cmd = args.app_path + "/bcftools concat -o " + args.out+"/somatic/" +  id + ".cell.gl.vcf.gz " +  tp + " -O z"
			print(cmd)
			output = os.system(cmd)

		joblst = []
		for id in chr_lst:
			joblst.append(id+">"+args.out)
		with Pool(processes=args.nthreads) as pool:
			result = pool.map(vcf2mat, joblst)
		error_check(all = chr_lst, output = result, step = "cellScan:vcf2mat")



	if args.step=="LDrefinement" or args.step=="all":
		logger.info("Run LD refinement ...")

		joblst = []
		for id in chr_lst:
			joblst.append(id+">"+args.out+">"+args.app_path)
		with Pool(processes=args.nthreads) as pool:
			result = pool.map(LDrefinement, joblst)
		error_check(all = chr_lst, output = result, step = "LDrefinement")



def error_check(all, output, step):
		job_fail = 0
		for id in all:
			if id not in output:
				logger.error("In "+ step + " step " + id + " failed!")
				job_fail = job_fail + 1

		if job_fail > 0:
			logger.error("Failed! See instructions above.")
			exit(1)



def preProcess(args):
	logger.info("Performing data preprocess before variant calling...")
	print_parameters_given(args)

	assert os.path.isfile(args.bamFile), "The bam file {} cannot be found!".format(args.bamFile)
	out = args.out
	os.system("mkdir -p " + out )
	os.system("mkdir -p " + out +  "/Bam")

	sample=[]
	with open(args.bamFile) as f_in:
			for line in f_in:
				record = line.strip().split(",")
				sample.append(record[0])
				#logger.debug("Checking sample {}".format(record[0]))
				assert len(record)==2, "Every line has to have exactly 2 comma-delimited columns! Line with sample name {} does not satisify this requiremnt!".format(record[0])
				assert os.path.isfile(record[1]), "Bam file {} cannot be found!".format(record[1])
				assert os.path.isfile(record[1]+".bai"), "Bam file {} has not been indexed!".format(record[1])
				#assert os.path.isabs(record[1]), "Please use absolute path for bam file {}!".format(record[1])

	para_lst = []
	with open(args.bamFile) as f_in:
			for line in f_in:
				record = line.strip().split(",")
				logger.debug("PreProcessing sample {}".format(record[0]))
				for chr in range(1, 23):
					para_single  =  dict(chr = "chr" + str(chr), 
						out = args.out, id = record[0], 
						bamFile = record[1], 
						max_mismatch = args.max_mismatch,
						samtools = samtools)
					para_lst.append(para_single)
	with Pool(processes=args.nthreads) as pool:
		result = pool.map(BamFilter, para_lst)
	# output the bam file list 

	# generate postProcess bam files 
	for chr in range(1, 23):
		bamlist = open(args.out + "/Bam/chr" +  str(chr) +  ".filter.bam.lst","w")
		for s in sample:
			bamlist.write(args.out+"/Bam/"+s+"_chr"+str(chr)+".filter.bam\n")
		bamlist.close()
		

def main():
	parser = argparse.ArgumentParser(
		description="""Monopogen: SNV calling from single cell sequencing
		""",
		epilog=
		"""Typical modules: preProcess, germline, somatic
		""",
		formatter_class=argparse.RawTextHelpFormatter)

	subparsers = parser.add_subparsers(title='Available subcommands', dest="subcommand")
	
	# every subcommand needs user config file
	common_parser = argparse.ArgumentParser(add_help=False)

	parser_preProcess = subparsers.add_parser('preProcess', parents=[common_parser],
		help='Preprocess of bam files including removing reads with high alignment mismatches',
		formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser_preProcess.add_argument('-b', '--bamFile', required=True,
								help="The bam file for the study sample, the bam file should be sorted. If there are multiple samples, each row with each sample") 
	parser_preProcess.add_argument('-o', '--out', required= False,
								help="The output director")
	parser_preProcess.add_argument('-a', '--app-path', required=True,
								help="The app library paths used in the tool")
	parser_preProcess.add_argument('-m', '--max-mismatch', required=False, type=int, default=3,
								help="The maximal alignment mismatch allowed in one reads for variant calling")
	parser_preProcess.add_argument('-t', '--nthreads', required=False, type=int, default=1,
								help="Number of threads used for SNVs calling")
	parser_preProcess.set_defaults(func=preProcess)


	parser_germline = subparsers.add_parser('germline', parents=[common_parser],
		help='Germline variant discovery, genotype calling from single cell data',
		formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser_germline.add_argument('-r', '--region', required= True, 
								help="The genome regions for variant calling")
	parser_germline.add_argument('-s', '--step', required= True, default="all",
								choices=['varScan', 'varImpute' , 'varPhasing', 'all'],
								help="Run germline variant calling step by step")
	parser_germline.add_argument('-o', '--out', required= False,
								help="The output director")
	parser_germline.add_argument('-g', '--reference', required= True, 
								help="The human genome reference used for alignment")
	parser_germline.add_argument('-p', '--imputation-panel', required= True, 
								help="The population-level variant panel for variant imputation refinement, such as 1000 Genome 3")
	parser_germline.add_argument('-m', '--max-softClipped', required=False, type=int, default=1,
								help="The maximal soft-clipped allowed in one reads for variant calling")
	parser_germline.add_argument('-a', '--app-path', required=True,
								help="The app library paths used in the tool")
	parser_germline.add_argument('-t', '--nthreads', required=False, type=int, default=1,
								help="Number of jobs used for SNVs calling")
	parser_germline.add_argument('-n', '--norun', required=False, default="FALSE", 
								choices=['TRUE','FALSE'],
								help="Generate the job scripts only. The jobs will not be run.")
	parser_germline.set_defaults(func=germline)

	parser_somatic = subparsers.add_parser('somatic', parents=[common_parser],
		help='Somatic variant calling from single cell sequencing',
		formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser_somatic.add_argument('-i', '--input-folder', required=True,
								help="The output folder from previous germline module")
	parser_somatic.add_argument('-r', '--region', required= True, 
								help="The chromosome IDs for variant calling. Each chromosomes in one row.")
	parser_somatic.add_argument('-l', '--barcode', required= True, 
								help="The csv file including cell barcode information")
	parser_somatic.add_argument('-k', '--keep', required= False, default=0.8,
								help="The proportion of reads kept for somatic calling. The cell will be sorted based on reads detected and cells with fewer reads will be removed.")
	parser_somatic.add_argument('-a', '--app-path', required=True,
								help="The app library paths used in the tool")
	parser_somatic.add_argument('-t', '--nthreads', required=False, type=int, default=22,
								help="Number of jobs used for SNV calling") 
	parser_somatic.add_argument('-w', '--winSize', required=False,  default="50MB",
								choices=['10MB','50MB'],
								help="Split the chromosome into small segments for cell level sequencing information collection. Setting 10MB will generate more jobs but be faster") 
	parser_somatic.add_argument('-s', '--step', required=True,
								choices=['featureInfo', 'cellScan' , 'LDrefinement', 'monovar', 'all'],
								help="Run germline variant calling step by step")
	parser_somatic.add_argument('-g', '--reference', required= True, 
								help="The human genome reference used for alignment")
	parser_somatic.set_defaults(func=somatic)

	args = parser.parse_args()

	if args.subcommand is None:
		# if no command is specified, print help and exit
		print("Please specify one subcommand! Exiting!")
		print("-"*80)
		parser.print_help()
		exit(1)

	# execute subcommand-specific function

	
	if args.subcommand == "somatic":
		args.out = args.input_folder

	global out, samtools, bcftools, bgzip, java, beagle 
	out = os.path.abspath(args.out)
	samtools  = os.path.abspath(args.app_path) + "/samtools" 
	bcftools = os.path.abspath(args.app_path) + "/bcftools"
	bgzip = os.path.abspath(args.app_path) + "/bgzip"
	java =  "java"
	beagle = os.path.abspath(args.app_path) + "/beagle.27Jul16.86a.jar"


	args.func(args)

	logger.info("Success! See instructions above.")



if __name__ == "__main__":
	main()
