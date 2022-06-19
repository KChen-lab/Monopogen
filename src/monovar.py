#!/usr/bin/env python
"""
The MIT License

Copyright (c) 2015
The University of Texas MD Anderson Cancer Center
Hamim Zafar and Ken Chen (kchen3@mdanderson.org)

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

"""

import sys
import re
import pysam
import time
import fileinput
import copy
import json
import _pickle as cPickle
import math
import copyreg as copy_reg
import types
import multiprocessing as mp
from functools import partial
from operator import add
from contextlib import closing
from utils import Utils_Functions
from alleles_prior import allele_prior
from Single_Cell_Ftrs_Pos import Single_Cell_Ftrs_Pos
from nu_prob_mat import Prob_matrix
from calc_variant_prob import Calc_Var_Prob
from genotype_prob_mat import Genotype_Prob_matrix
from nu_genotype_single_cell import Single_cell_genotype_records
from mp_genotype import MP_single_cell_genotype
from hzvcf import VCFRecord, VCFDocument, VRecord, VCF

# Required for using multiprocessing


def _pickle_method(m):
    if m.im_self is None:
        return getattr, (m.im_class, m.im_func.func_name)
    else:
        return getattr, (m.im_self, m.im_func.func_name)

copy_reg.pickle(types.MethodType, _pickle_method)

U = Utils_Functions()
M = MP_single_cell_genotype()

# Default values
pe = 0.002
pad = 0.2
thr = 0.05
m_thread = 1
CF_flag = 1

input_args = {}

# Process the inputs
argc = len(sys.argv)
i = 1
while (i < argc):
    if (sys.argv[i] == '-n'):
        n_cells = int(sys.argv[i + 1])      # Number of input bam files
    elif (sys.argv[i] == '-p'):
        pe = float(sys.argv[i + 1])    	  # probability of error
    elif (sys.argv[i] == '-d'):
        pd = float(sys.argv[i + 1])    	  # probability of deamination error
    elif (sys.argv[i] == '-a'):
        pad = float(sys.argv[i + 1])  	  # Probability of ADO
    elif (sys.argv[i] == '-f'):
        ref_file = sys.argv[i + 1]          # Reference Genome File
        input_args['-f'] = 'Provided'
    elif (sys.argv[i] == '-b'):
        bam_file_list = sys.argv[i + 1]	  # File containing list of bam files
        input_args['-b'] = 'Provided'
    elif (sys.argv[i] == '-o'):
        outfile = sys.argv[i + 1]      	  # Output File
        input_args['-o'] = 'Provided'
    elif (sys.argv[i] == '-t'):
        # threshold to use for calling variant
        thr = float(sys.argv[i + 1])
    elif (sys.argv[i] == '-c'):
        CF_flag = int(sys.argv[i + 1])      # Flag for using Consensus Filter
    elif (sys.argv[i] == '-m'):
        # Number of threads to use in multiprocessing
        m_thread = int(sys.argv[i + 1])
    i = i + 2

 
try:
    b = input_args['-f']
except KeyError:
    print ("Error: Reference genome file not provided. Use '-f' for reference genome file.\n")
    exit(3)

try:
    b = input_args['-b']
except KeyError:
    print ("Error: List of Bam files not provided. Use '-b' for list of Bam files.\n")
    exit(3)

try:
    b = input_args['-o']
except KeyError:
    print ("Error: Output file not provided. Use '-o' for Output file.\n")
    exit(3)
try:
    assert CF_flag <= 1
except AssertionError:
    print ("CF_flag can have value 0 or 1. Use '-c' with proper value.\n")
    exit(3)

# Obtain the RG IDs from the bam files
bam_id_list = []
f_bam_list = open(bam_file_list)
for filename in f_bam_list:
    filename = filename.replace('\n', '')
    print(filename)
    bam_id = U.Get_BAM_RG(filename)
    bam_id_list.append(bam_id)

n_cells = len(bam_id_list)

# Initialize the pool of multiprocessing
pool = mp.Pool(processes=m_thread)

# Constants to be used later
cell_no_threshold = n_cells / 2
# no of possible alternate alleles {0, 1, 2, ..., 2m}
max_allele_cnt = 2 * n_cells + 1
theta = 0.001  # Heterozygosity rate
Base_dict = {0: 'A', 1: 'T', 2: 'G', 3: 'C'}

# Factorial values and nCr values
# Table for all the required factorials
factorial_list = U.Create_Factorial_List(max_allele_cnt)
# Table for all the required nCr
nCr_matrix = U.Create_nCr_mat(max_allele_cnt, factorial_list)

# Dictionary for holding all the priors for different values of n
prior_variant_dict = {}
for i in range(n_cells + 1):
    prior_variant_dict[i] = U.calc_prior(theta, i, 1)

# Open VCF file and print header
f_vcf = open(outfile, 'w')
vcf = VCFDocument(f_vcf)
vcf.populate_fields(bam_id_list)
vcf.populate_reference(ref_file)
vcf.print_header()

# List of all single_cell_ftr_pos object
All_single_cell_ftrs_list = n_cells * [None]
# Global list for storing which cell contains read support
read_flag_row = n_cells * [None]
# Global list for storing which cell has alternate allele support
alt_allele_flag_row = n_cells * [None]

for line in sys.stdin:
    line = line.replace('\n', '')
    row = line.split('\t')
    contig = row[0]
    pos = int(row[1])
    refBase = U.refineBase(row[2])
    total_depth = 0
    total_ref_depth = 0
    for i in range(1, n_cells + 1):
        curr_cell_pos_ftrs = Single_Cell_Ftrs_Pos(
            refBase, row[3 * i: 3 * i + 3])
        total_depth += curr_cell_pos_ftrs.depth
        total_ref_depth += curr_cell_pos_ftrs.refDepth
        All_single_cell_ftrs_list[i - 1] = curr_cell_pos_ftrs

    Alt_count = total_depth - total_ref_depth
    if total_depth == 0:
        Alt_freq = 0
    else:
        Alt_freq = float(Alt_count) / total_depth

    # No reads supporting alternate allele, so no operations needed
    if (total_ref_depth == total_depth):
        continue

    # Cases that are to be prefiltered
    elif ((total_depth > 30) & ((Alt_count <= 2) | (Alt_freq <= 0.001))):
        continue

    elif (Alt_freq <= 0.01):
        continue

    # Bad reference, so filtered
    elif (refBase not in ['A', 'T', 'G', 'C']):
        continue

    elif (total_depth <= 10):
        continue

#	elif ((total_depth > 1000) & (Alt_freq <= 0.003)):
#		continue

    else:
        # List for storing the sngl_cell_objs that have read support, will be
        # further used in the model
        read_supported_cell_list = []
        # Gloabal list for storing the alternate allele counts
        total_alt_allele_count = [0, 0, 0, 0]
        # Global list for storing the indices of the cells having read support
        cell_index_list = []
        info_list = ['GT:AD:DP:GQ:PL']
        max_depth = 10000

        # Traverse through all the sngl_cell_ftr_obj and if has read support
        # further calculate the other quantities
        c = 1
        for j in range(n_cells):
            sngl_cell_ftr_obj = All_single_cell_ftrs_list[j]
            read_flag = U.checkReadPresence(sngl_cell_ftr_obj)
            if read_flag == 1:
                alt_allele_flag = U.CheckAltAllele(sngl_cell_ftr_obj)
                sngl_cell_ftr_obj.Get_base_call_string_nd_quals(refBase)
                sngl_cell_ftr_obj.Get_Alt_Allele_Count()   # Get the alt allele counts
                sngl_cell_ftr_obj.Set_Cell_Index(
                    c)        # Save the cell index
                if alt_allele_flag == 1:
                    # Update the list of total_alt_allele_count
                    total_alt_allele_count = U.update_alt_count(
                        total_alt_allele_count, sngl_cell_ftr_obj)
                # Populate the list of read supported cells
                read_supported_cell_list.append(sngl_cell_ftr_obj)
                # Populate cell_index_list
                cell_index_list.append(c)
            else:
                alt_allele_flag = 0
            read_flag_row[j] = read_flag
            alt_allele_flag_row[j] = alt_allele_flag
            c = c + 1

        # Number of cells having read support
        read_smpl_count = sum(read_flag_row)
        # Number of cells having alternate allele support
        alt_smpl_count = sum(alt_allele_flag_row)
        Alt_count = max(total_alt_allele_count)    # Update the Alt_count value
        # Get the altBase
        if Alt_count > 0:
            altBase = Base_dict[total_alt_allele_count.index(Alt_count)]
        else:
            altBase = ''

        if (altBase == ''):
            continue

        # Calculate prior_allele_mat
        prior_allele_mat = U.Get_prior_allele_mat(
            read_smpl_count, alt_smpl_count, cell_no_threshold, total_depth, Alt_freq, pe)

        # Operations on the single cells with read support
        # Number of cells with read support
        read_supported_n_cells = len(read_supported_cell_list)
        for j in range(read_supported_n_cells):
            read_supported_cell_list[j].Store_Addl_Info(
                refBase, altBase, Alt_freq, prior_allele_mat)

        # Get prior_variant_number distribution
        prior_variant_number = prior_variant_dict[read_supported_n_cells]

        if (read_supported_n_cells == 0):
            continue

        # Obtain the value of probability of SNV
        Calc_var_prob_obj = Calc_Var_Prob(
            read_supported_cell_list, prior_allele_mat, read_supported_n_cells)
        (zero_variant_prob, denominator) = Calc_var_prob_obj.Calc_Zero_Var_Prob(
            n_cells, max_depth, nCr_matrix, pad, prior_variant_number)

        # Probability of SNV passes the threshold
        if zero_variant_prob <= thr:

            (max_prob_ratio, max_prob_allele_count) = U.find_max_prob_ratio(
                Calc_var_prob_obj.matrix, Calc_var_prob_obj.matrix_shape)
            (oddsRatio, pval) = U.calc_strand_bias(
                read_supported_cell_list, Alt_count)

            if (total_ref_depth > 0):
                (baseQranksum, baseQranksum_p_val) = U.Calc_Base_Q_Rank_Sum(
                    read_supported_cell_list, refBase, altBase)
            else:
                baseQranksum = 0.0
            genotype_dict = {0: '0/0', 1: '0/1', 2: '1/1'}
            barcode = '<'
            func = partial(M.get_info_string, read_supported_cell_list, prior_allele_mat,
                           n_cells, nCr_matrix, prior_variant_number, denominator, genotype_dict)
            output = pool.map(func, range(read_supported_n_cells))
            read_supported_info_list = [p[0] for p in output]
            read_supported_barcodes = [p[1] for p in output]
            for j in range(n_cells):
                if (All_single_cell_ftrs_list[j].depth == 0):
                    info_list.append('./.')
                    barcode += 'X'
                else:
                    info_list.append(read_supported_info_list[0])
                    del read_supported_info_list[0]
                    barcode += read_supported_barcodes[0]
                    del read_supported_barcodes[0]
            barcode += '>'
            if zero_variant_prob == 0:
                Qual = 323
            else:
                Qual = 0 - math.log10(abs(zero_variant_prob))

            (AC, AF, AN) = U.Calc_chr_count(barcode)
            QD = U.Calc_qual_depth(barcode, All_single_cell_ftrs_list, Qual)
            PSARR = U.Calc_Per_Smpl_Alt_Ref_Ratio(
                total_ref_depth, Alt_count, read_smpl_count, alt_smpl_count)
            vcf_record = VRecord(contig, pos)
            info_record = [str(AC), "%.2f" % AF, str(AN), "%.2f" % baseQranksum, str(
                total_depth), str(QD), "%.2f" % oddsRatio, "%.2f" % max_prob_ratio, "%.2f" % PSARR]
            info_names = [i[0] for i in vcf.info_fields]
            info_string = ';'.join(
                '%s=%s' % t for t in zip(info_names, info_record))

            if CF_flag == 1:
                filter_flag = U.Consensus_Filter(barcode)
                if filter_flag == 1:
                    vcf_record.get6fields(
                        refBase, altBase, '.', Qual, 'PASS', info_string)
                else:
                    vcf_record.get6fields(
                        refBase, altBase, '.', Qual, '.', info_string)
                vcf_record.format_vcf(info_list)
                vcf_record.get_passcode(barcode)
                vcf.print_my_record(vcf_record)

            else:
                vcf_record.get6fields(
                    refBase, altBase, '.', Qual, '.', info_string)
                vcf_record.format_vcf(info_list)
                vcf_record.get_passcode(barcode)
                vcf.print_my_record(vcf_record)
