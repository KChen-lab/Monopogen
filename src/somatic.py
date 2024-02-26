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
from bamProcess import * 
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
	print(cmd)
	os.system(cmd)
	#process = subprocess.run(cmd, shell=True)


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
			phased_vcf = args.out + "/germline/" + jobid + ".phased.vcf.gz" 
			assert os.path.isfile(bam_filter), "The bam list file {} cannot be found! Please run germline module".format(bam_filter)
			assert os.path.isfile(gl_vcf), "The *.gl.vcf.gz file {} cannot be found! Please run germline module".format(gl_vcf)
			assert os.path.isfile(phased_vcf), "The *.phased.vcf.gz file {} cannot be found! Please run germline module".format(phased_vcf)
			assert os.path.isfile(args.barcode), "The cell barcode file {} cannot be found!".format(args.barcode)


def getInfo_robust(rec, info):

	info_dt =  rec.info.get(info)
	if type(info_dt) is not type(None):
		if not isinstance(info_dt,float):
			info_dt = info_dt[0]
		info_dt = round(info_dt,2)
	return(info_dt)


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

	#os.system("cat " +  out + "/somatic/" + chr + "*.gl.vcf.filter.hc.bed > " + out + "/somatic/" + chr + ".bed")
	chr_bed = out + "/somatic/" + chr + ".gl.vcf.filter.hc.bed"
	cmd1 = samtools + " view " + " -b  -L " + chr_bed + " " + inbam + " -o " + outbam
	with open(out+"/Script/bamExtract_" + chr + ".sh","w") as f_out:
		f_out.write(cmd1 + "\n")
		f_out.write(samtools + " index " +  outbam + "\n")
	cmd="bash " + out+"/Script/bamExtract_" + chr + ".sh"
	#print(cmd)
	os.system(cmd)

def featureInfo(para):

	para_lst = para.strip().split(">")
	region = para_lst[0]
	out = para_lst[1]
	app = para_lst[2]
	pysam.tabix_index(out + "/germline/" +  region + ".phased.vcf.gz", preset="vcf", force=True)
	pysam.tabix_index(out + "/germline/" +  region + ".gl.vcf.gz", preset="vcf", force=True)
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
		info_VDB = getInfo_robust(rec,'VDB')
		info_RPB = getInfo_robust(rec,'RPB')
		info_MQB = getInfo_robust(rec,'MQB')  
		info_BQB = getInfo_robust(rec,'BQB') 
		info_MQSB =getInfo_robust(rec,'MQSB')
		info_SGB = getInfo_robust(rec,'SGB')
		info_MQ0F =getInfo_robust(rec,'MQ0F')
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
	gl_vcf_filter_dp4.close()
	gl_vcf_filter_txt.close() 
	gl_vcf_filter_bed.close() 
	gl_vcf_dp4.close()

	
	bamExtract(para)
	return(region)


def getBamName(chr, out):
	#chr = region.split(":")[0]
	infile = out + "/Bam/" + chr + ".filter.bam.lst"
	with open(infile) as f_in:
		for line in f_in:
			bamName = line.strip()
			return(bamName)



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






def robust_get_tag(read, tag_name):  
	try:  
		return read.get_tag(tag_name)
	except KeyError:
		return "NotFound"

def rev_compl(st):
    nn = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(nn[n] for n in reversed(st))



def bam2mat(para):

	para_lst = para.strip().split(">")
	#print(para_lst)

	bam2gz(para);

	#### gz to matrix ### 
	mat_infile = para_lst[1] + "/somatic/" + para_lst[0] + ".cell_snv.mat.gz"
	snv_infile = para_lst[1] + "/somatic/" + para_lst[0] + ".cell_snv.snvID.csv"
	cell_infile = para_lst[1] + "/somatic/" + para_lst[0] + ".cell_snv.cellID.csv"
	meta_info = para_lst[1] + "/somatic/" + para_lst[0] + ".gl.vcf.filter.DP4"

	cell_clst = pd.read_csv(cell_infile)
	snv_clst =  pd.read_csv(snv_infile)  
	mat = pd.read_csv(mat_infile, sep="\t", header=None) 
	mat_out = gzip.open(para_lst[1] + "/somatic/" +  para_lst[0] + ".gl.filter.hc.cell.mat.gz","wt")


	mat.columns = ['snvIndex','cellIndex','allele'] 
	mat=mat.groupby(by=['snvIndex','cellIndex'], as_index=False).first()
	mat=mat.pivot(index='snvIndex', columns='cellIndex', values='allele')
	meta_info = pd.read_csv(meta_info, sep="\t", header=None) 
	meta_info.columns=["chr","pos","ref_allele","alt_allele","Dep","dep1","dep2","dep3","dep4","genotype","QS","VDB","RPB","MQB","BQB","MQSB","SGB","MQ0F"]

	overlapSNV_index=set(snv_clst["index"]).intersection(set(mat.index))

	cell_clst.loc[list(mat.columns)]['cell'].to_csv(para_lst[1] + "/somatic/" + para_lst[0] + ".cell_snv.cellID.filter.csv")

	## meta_info_id is similar with snv_clst
	n_cell = mat.shape[1]
	for index in overlapSNV_index:
		info = meta_info.iloc[index].astype(str)
		geno = info["genotype"]
		info = "\t".join(info)
		
		flag = 0
		if geno=="0|1":
			phase_info = ["0|0"]*n_cell
			tp = mat.loc[index]
			pos = tp.dropna().index
			for value in pos:
					allele = tp.loc[value]
					index_vec = tp.index.get_loc(value)
					if allele==1:
						phase_info[index_vec]="0|1"
					elif allele==0:
						phase_info[index_vec]="1|0"
			flag = 1
		elif geno=="1|0":
			phase_info = ["0|0"]*n_cell
			tp = mat.loc[index]
			pos = tp.dropna().index
			for value in pos:
					allele = tp.loc[value]
					index_vec = tp.index.get_loc(value)
					if allele==1:
						phase_info[index_vec]="1|0"
					elif allele==0:
						phase_info[index_vec]="0|1"
			flag = 1 
		elif geno=="nan":
			phase_info = ["0/0"]*n_cell
			tp = mat.loc[index]
			pos = tp.dropna().index
			for value in pos:
					allele = tp.loc[value]
					index_vec = tp.index.get_loc(value)
					if allele==1:
						phase_info[index_vec]="0/1"
					elif allele==0:
						phase_info[index_vec]="1/0"
			flag = 1 
		if flag:
			mat_out.write(info+"\t" +'\t'.join(phase_info))
			mat_out.write('\n')
	return(para_lst[0])




def bam2gz(para):
			
		para_lst = para.strip().split(">")

		#assert os.path.isfile(bam_filter), "*.fiter.targeted.bam file {} cannot be found!".format(bam_filter)
		#in_fasta = "/rsrch3/scratch/bcb/jdou1/scAncestry/ref/fasta/genome.fa"
		#in_vcf = "/rsrch3/scratch/bcb/jdou1/bam_process/chr1.gp.vcf"
		#in_bam = "/rsrch3/scratch/bcb/jdou1/bam_process/chr1.filter.targeted.bam"

		# joblst.append(id+">"+args.out+">"+args.reference)
		
		in_snv = para_lst[1] + "/somatic/" + para_lst[0] + ".gl.vcf.filter.DP4"
		in_bam = para_lst[1] + "/Bam/" + para_lst[0] + ".filter.targeted.bam"
		in_fasta = para_lst[2]
		in_cell_barcode = para_lst[3]

		# read cell barcode file 

		cell_clst = pd.read_csv(in_cell_barcode)   
		df = pd.DataFrame(cell_clst, columns= ['cell','id'])
		df = df.sort_values(by=['id'], ascending=False)
		df['index']=range(len(df.index))
		cell_set = list(df["cell"])


		ref_fa = pysam.FastaFile(in_fasta)
		snv_info = {}
		tp = list()
		index = 0 
		motif_len = 3
		snv_tol = 0 

		# index 0-based 
		with open(in_snv,'r') as fp:
			for line in fp:
				line = line.strip()
				line = line.split("\t")
				chrom = line[0]
				pos = int(line[1])
				ref = line[2] 
				alts = line[3]
				seq = ref_fa.fetch(chrom, pos-4, pos+3)
				seq_rev_compl = rev_compl(seq)
				id = str(chrom)+":"+str(pos) + ":" + ref + ":" + alts[0]
				mydict = dict(id=index, chr = chrom, pos = pos, motif_pos= seq, ref_allele = ref, alt_allele=alts[0])			
				snv_info[index]=mydict  
				index = index + 1
				snv_tol = snv_tol + 1
				tp.append(id)

		

		#print(snv_info_out)
		## output cell snv meta information 
		df.to_csv(para_lst[1] + "/somatic/" + para_lst[0] + ".cell_snv.cellID.csv", index = False) 
		snv_info_out = pd.DataFrame()
		snv_info_out['snvID'] = tp 
		snv_info_out['index'] = range(len(tp))
		snv_info_out.to_csv(para_lst[1] + "/somatic/" + para_lst[0] + ".cell_snv.snvID.csv",index = False) 
		
		read_tol = 0 
		read_cover_tol = 0 
		read_wild_tol = 0
		read_mutated_tol = 0 
		read_noAllele_tol = 0 
		overlap = 0 
		index = 0
		lower_index = 0

		infile = pysam.AlignmentFile(in_bam,"rb")
		fp = gzip.open(para_lst[1] + "/somatic/" + para_lst[0] + ".cell_snv.mat.gz", "wt")
		# fp = open(para_lst[1] + "/somatic/" + para_lst[0] + ".cell_snv.mat", "wt")
		for s in infile:

			cell_barcode  = robust_get_tag(s,"CB")
			if cell_barcode=="NotFound":
				### no cell barcode file detected.
				### check whether the barcode added in the query name 	
				cell_barcode = s.query_name.strip().split("_")[1]

			if cell_barcode in cell_set:
				cell_barcode_index = cell_set.index(cell_barcode)
			else:
				### cell barcode is not detected, should move to next read 
				continue

			read_name = s.query_name 
			align_chr = s.reference_name
			align_start = s.reference_start
			align_seq = s.query_sequence
			mystart = str(snv_info[lower_index]["pos"])


			### flags used to see whether SNVs covered in different scenarios
			read_tol = read_tol + 1 
			read_cover = 0
			read_wild = 0 
			read_mutated = 0
			read_noAllele = 0
			read_len = len(align_seq)


			### get the SNVs locating in [align_start, align_start + read_len] ### 
			if read_tol%1000000 == 0:
				print("scanning read " + str(read_tol)) 

			lock = 0 
			snv_cover = ""


			for i in range(lower_index,snv_tol-1,1):

				snv_pos = snv_info[i]["pos"]
				wild_allele = snv_info[i]["ref_allele"]
				alt_allele = snv_info[i]["alt_allele"]

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
							fp.write(str(i)+"\t"+str(cell_barcode_index)+"\t0\n")
						elif re.search(motif_neg, align_seq):
							read_mutated = read_mutated + 1
							fp.write(str(i)+"\t"+str(cell_barcode_index)+"\t1\n")
						else:
							delta = snv_pos - align_start	
							if read_len - delta < motif_len:
								motif_pos_part = motif_pos[0:motif_len+1]
								motif_neg_part = motif_neg[0:motif_len+1]
								seq_part = align_seq[read_len-2*motif_len-1:read_len]
								
								if re.search(motif_pos_part, seq_part):
									read_wild = read_wild + 1 	
									fp.write(str(i)+"\t"+str(cell_barcode_index)+"\t0\n")		
								elif re.search(motif_neg_part, seq_part):	
									read_mutated = read_mutated + 1
									fp.write(str(i)+"\t"+str(cell_barcode_index)+"\t1\n")

							elif delta <= motif_len: 
								motif_pos_part = motif_pos[motif_len:len(motif_pos)]
								motif_neg_part = motif_neg[motif_len:len(motif_neg)]
								seq_part = align_seq[0:2*motif_len+1]

								if re.search(motif_pos_part, seq_part):
									read_wild = read_wild + 1 
									fp.write(str(i)+"\t"+str(cell_barcode_index)+"\t0\n")	
								elif re.search(motif_neg_part, seq_part):
									read_mutated = read_mutated + 1
									fp.write(str(i)+"\t"+str(cell_barcode_index)+"\t1\n")
							else:
								#fp.write(str(snv_info[i])+"\n")
								#fp.write(str(align_chr) + ":" + str(align_start) + ":" + align_seq+"\n") 
								read_noAllele = read_noAllele + 1 
					else:
						break 
						
			if read_cover > 0: 
				read_cover_tol = read_cover_tol + 1 
			if read_wild > 0: 
				read_wild_tol = read_wild_tol + 1 
			if read_mutated > 0 and read_wild>0: 
				overlap = overlap + 1

			if read_mutated > 0: 
				read_mutated_tol = read_mutated_tol + 1 
			if read_noAllele > 0: 
				read_noAllele_tol = read_noAllele_tol + 1


		infile.close()
		fp.close()	
		print(str(read_tol) + ":" + str(read_cover_tol) + ":" + str(read_wild_tol) + ":" + str(read_mutated_tol) + ":" + str(overlap) + ":" + str(read_noAllele_tol))
		return(para_lst[0])