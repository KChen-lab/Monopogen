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



def withSNVs(invcf, path):
	print(invcf)
	#pysam.tabix_index(invcf, preset="vcf", force=TRUE)
	print(path)
	os.system(path + "/tabix -p vcf " + invcf)
	vcf = VariantFile(invcf) 
	cnt = 0
	for rec in vcf.fetch():
		cnt = cnt + 1 
	return(cnt)

def runCMD(cmd):

	os.system(cmd)
	#process = subprocess.run(cmd, shell=True, stdout=open(args.logfile, 'w'), stderr=open(args.logfile,'w'))


def robust_get_tag(read, tag_name):  
	try:  
		return read.get_tag(tag_name)
	except KeyError:
		return "NotFound"


def validate_user_setting_somatic(args):

	assert os.path.isdir(args.out), "The germline output folder {} cannot be found! Please run germline module.".format(args.out)
	assert os.path.isfile(args.region), "The region file {} cannot be found! Please set the genomic regions for somatic detection".format(args.region)
	
	with open(args.region) as f_in:
		for line in f_in:
			record = line.strip().split(",")
			if(len(record)==1):
				jobid = record[0]
			if(len(record)==3):
				jobid = record[0] + ":" + record[1] + "-" + record[2]
				logger.error("Only the whole chromsome calling is allowed!")
				exit(1)
			bam_filter = args.out + "/Bam/" +  record[0] +  ".filter.bam.lst"
			gl_vcf = args.out + "/germline/" + jobid + ".gl.vcf.gz"
			gp_vcf = args.out + "/germline/" + jobid + ".germline.vcf.gz" 
			phased_vcf = args.out + "/germline/" + jobid + ".phased.vcf.gz" 
			assert os.path.isfile(bam_filter), "The bam list file {} cannot be found! Please run germline module".format(bam_filter)
			assert os.path.isfile(gl_vcf), "The *.gl.vcf.gz file {} cannot be found! Please run germline module".format(gl_vcf)
			#if withSNVs(gl_vcf, args.app_path)==0:
			#	print("The *.gl.vcf.gz file {} cannot be found! Maybe there are no markers detected in ".format(gl_vcf))
			#if withSNVs(phased_vcf, args.app_path)==0:
			#	print("The *.phased.vcf.gz file {} cannot be found! Maybe there are no markers detected in this region?".format(phased_vcf))
			assert os.path.isfile(args.barcode), "The cell barcode file {} cannot be found!".format(args.barcode)


def getInfo_robust(rec, info):

	info_dt =  rec.info.get(info)
	if type(info_dt) is not type(None):
		if not isinstance(info_dt,float):
			info_dt = info_dt[0]
		info_dt = round(info_dt,2)
	return(info_dt)


def featureInfo(para):
	para_lst = para.strip().split(">")
	out = para_lst[1]
	region = para_lst[0]

	vcf_in = VariantFile(out + "/germline/" +  region + ".phased.vcf.gz") 
	info_GT = {}
	for rec in vcf_in.fetch():
		GT = [value['GT'] for value in rec.samples.values()][0]
		id = str(rec.chrom)+":"+str(rec.pos) + ":" + rec.ref + ":" + rec.alts[0]
		info_GT[id] = GT

	gl_vcf_in = VariantFile(out + "/germline/" +  region + ".gl.vcf.gz") 
	gl_vcf_dp4 = open(out + "/somatic/" +  region + ".gl.vcf.DP4","w")
	gl_vcf_filter_dp4 = open(out + "/somatic/" +  region + ".gl.vcf.filter.DP4","w")
	gl_vcf_filter_bed = open(out + "/somatic/" +  region + ".gl.vcf.filter.hc.bed","w")
	gl_vcf_filter_txt = open(out + "/somatic/" +  region + ".gl.vcf.filter.hc.pos","w")

	depth_filter_novelSNV = 10
	for rec in gl_vcf_in.fetch():
		info_I16 = rec.info.get('I16')
		info_QS =  getInfo_robust(rec,'QS')
		info_VDB =  getInfo_robust(rec,'VDB')
		info_RPB = getInfo_robust(rec,'RPB')
		info_MQB = getInfo_robust(rec,'MQB')  
		info_BQB = getInfo_robust(rec,'BQB') 
		info_MQSB = getInfo_robust(rec,'MQSB')
		info_SGB = getInfo_robust(rec,'SGB')
		info_MQ0F = getInfo_robust(rec,'MQ0F')
		id = str(rec.chrom)+":"+str(rec.pos) + ":" + rec.ref + ":" + rec.alts[0]
		gt_info_var = "NA"
		if (id in info_GT):
			gt_info_var = str(info_GT[id][0]) + "|" + str(info_GT[id][1])
		a = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(rec.chrom, rec.pos, 
			rec.ref, rec.alts[0], rec.info['DP'], 
			info_I16[0], info_I16[1], info_I16[2], info_I16[3], gt_info_var, info_VDB, info_QS, info_RPB,
			info_MQB, info_BQB, info_MQSB, info_SGB, info_MQ0F)
		gl_vcf_dp4.write(a)

		# include all variants from germline 
		if (info_I16[0] + info_I16[1]>=4 and info_I16[2] + info_I16[3]>=4 or (id in info_GT)):
			b = "{}\t{}\t{}\n".format(rec.chrom, rec.pos-1, rec.pos)
			gl_vcf_filter_bed.write(b)
			b = "{}\t{}\n".format(rec.chrom, rec.pos)
			gl_vcf_filter_txt.write(b)
			gl_vcf_filter_dp4.write(a)

def getBamName(chr, out):
	#chr = region.split(":")[0]
	infile = out + "/Bam/" + chr + ".filter.bam.lst"
	with open(infile) as f_in:
		for line in f_in:
			bamName = line.strip()
			return(bamName)

def bamExtract(para):
	
	para_lst = para.strip().split(">")
	
	chr = para_lst[0]
	out = para_lst[1]
	app_path = para_lst[2]
	samtools = os.path.abspath(app_path) + "/samtools" 
	out = os.path.abspath(out)
	inbam = getBamName(chr, out)
	outbam =  out + "/Bam/" + chr + ".filter.targeted.bam"
	# remember to update bed file   ls -lrt

	os.system("cat " +  out + "/somatic/" + chr + "*.gl.vcf.filter.hc.bed > " + out + "/somatic/" + chr + ".bed")
	cmd1 = samtools + " view " + " -b  -L " + out + "/somatic/" + chr +".bed " + inbam + " -o " + outbam
	with open(out+"/Script/bamExtract_" + chr + ".sh","w") as f_out:
		f_out.write(cmd1 + "\n")
		f_out.write(samtools + " index " +  outbam + "\n")
	cmd="bash " + out+"/Script/bamExtract_" + chr + ".sh"
	runCMD(cmd)


def bamSplit(para):
	para_lst = para.strip().split(":")
	chr = para_lst[0]
	cell = para_lst[1]
	out = para_lst[2]
	app_path = para_lst[3]
	bam_filter =  out + "/Bam/" + chr + ".filter.targeted.bam"
	#assert os.path.isfile(bam_filter), "*.fiter.targeted.bam file {} cannot be found!".format(bam_filter)
	samtools = app_path + "/samtools" 	
	infile = pysam.AlignmentFile(bam_filter,"rb")
	# Note to change the read groups 
	tp =infile.header.to_dict()
	if len(tp['RG'])>1:
		tp['RG']= [tp['RG'][0]]
	tp['RG'][0]['SM'] = cell
	tp['RG'][0]['ID'] = cell
	cnt = 0
	outfile =  pysam.AlignmentFile( out + "/Bam/split_bam/" + chr + "_" + cell  + ".bam", "wb", header=tp)
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
	cmd=samtools + " index " + out + "/Bam/split_bam/" + chr + "_" + cell + ".bam"
	#runCMD(cmd,args)
	os.system(samtools + " index " + out + "/Bam/split_bam/" + chr + "_" + cell + ".bam")

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
	bam_filter = out + "/Bam/split_bam/cell" + chr + ".bam.lst"
	cmd1 = samtools + " mpileup -b " + bam_filter + " -f "  + reference + " -r " +  jobid + " -q 20 -Q 20 -t DP4 -d 10000 -v "
	cmd1 = cmd1 + " | " + bcftools + " view " + " | "  + bcftools  + " norm -m-both -f " + reference 
	cmd1 = cmd1 + " | grep -v \"<X>\" | grep -v INDEL |" + bgzip +   " -c > " + out + "/somatic/" +  jobid + ".cell.gl.vcf.gz" 		
	print(cmd1)
	output = os.system(cmd1)

	# delete the bam files once snv calling was finished in specfic regions
	f = open(bam_filter, "r")
	for x in f:
		x = x.strip()
		os.system("rm " + x)
		os.system("rm " + x + ".bai")		
	f.close()

	if output == 0:
		return(jobid)






def less1(num):
	if num>1:
		num = 1 
	return num

def vcf2mat(para):
	para_lst = para.strip().split(">")
	region = para_lst[0]
	out = para_lst[1]
	vcf_in = pysam.VariantFile(out + "/somatic/" +  region + ".cell.gl.vcf.gz") 
	mat_out = gzip.open(out + "/somatic/" +  region + ".gl.filter.hc.cell.mat.gz","wt")
	meta_info = {}
	phase_info = {}
	with open(out + "/somatic/" +  region + ".gl.vcf.filter.DP4",'r') as fp:
		for line in fp:
			line = line.strip()
			data = line.split("\t")
			id=str(data[0])+":"+str(data[1])+":"+data[2] + ":" + data[3]
			meta_info[id] = line
			phase_info[id] = data[9]
	n = len(list((vcf_in.header.samples)))


	with open(out + "/somatic/" +  region + ".cell.txt",'wt') as fp:
		for line in list(vcf_in.header.samples):
			fp.write(line + "\n")
	fp.close()

	var_lst = {}
	cnt = 0
	for rec in vcf_in.fetch():
		#print(rec.samples.values())
		allele =  rec.ref + ":" + rec.alts[0]
		id = str(rec.chrom)+":"+str(rec.pos)+ ":" + allele
		var_lst[id] = 1 
		a = [None]*(n+1)
		if (id in meta_info):
			a[0] = meta_info[id]
			i = 0
			cnt = cnt + 1
			if phase_info[id]=="0|1": 
				for value in rec.samples.values():
					ref = value['DP4'][0] +  value['DP4'][1]
					alt = value['DP4'][2] +  value['DP4'][3]
					ref = less1(ref)
					alt = less1(alt)
					i = i + 1
					a[i] = str(ref)+"|"+str(alt)
				mat_out.write('\t'.join(a))
				mat_out.write('\n')
			if phase_info[id]=="1|0": 
				for value in rec.samples.values():
					ref = value['DP4'][0] +  value['DP4'][1]
					alt = value['DP4'][2] +  value['DP4'][3]
					ref = less1(ref)
					alt = less1(alt)
					i = i + 1
					a[i] = str(alt)+"|"+str(ref)
				mat_out.write('\t'.join(a))
				mat_out.write('\n')
			if phase_info[id]=="NA":
				for value in rec.samples.values():
					ref = value['DP4'][0] +  value['DP4'][1]
					alt = value['DP4'][2] +  value['DP4'][3]
					ref = less1(ref)
					alt = less1(alt)
					i = i + 1
					a[i] = str(ref)+"/"+str(alt)
				mat_out.write('\t'.join(a))
				mat_out.write('\n')

	mat_out.close()
	vcf_in.close()
	return(region)
    
def LDrefinement(para):
	para_lst = para.strip().split(">")
	region = para_lst[0]
	out = para_lst[1]
	app_path = para_lst[2]
	outdir = out+"/somatic/"
	cellfile = out+"/somatic/"+region+".gl.filter.hc.cell.mat.gz"
	cmd = "Rscript " + app_path + "/../src/LDrefinement.R  " + cellfile + " " + outdir + " " + region 
	print(cmd)
	output = os.system(cmd)
	if output == 0:
		return(region)

