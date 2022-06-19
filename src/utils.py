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

import math
import heapq
import copy
import re
import sys
import numpy as np
import pysam
from operator import add
from scipy import stats
import fileinput
from contextlib import closing
from base_q_ascii import Ascii_Table
from alleles_prior import allele_prior

base_q_tbl = Ascii_Table()


class Utils_Functions:

    def nCr(self, n, r):
        f = math.factorial
        return f(n) / f(r) / f(n - r)

    def find_second_smallest(self, l):
        return heapq.nsmallest(2, l)[-1]

    def GetPosition(self, row):
        """Get the position in the chromosome from the line of the file
        """
        position = int(row[1])
        return position

    def GetRefBase(self, row):
        """ Get the reference base at that position from the line of the file
        """
        refbase = row[2]
        return refbase

    def single_cell_ftrs(self, row, cell_index):
        """ Extract the features of the cell indexed with cell_index
        """
        newrow = row[cell_index + 2]
        newrow = newrow.split('\t')
        cell_ftr = [int(i) for i in newrow[0:6]]
        read_bases = newrow[6]
        read_bases = read_bases.replace(' ', '')
        newrow[7] = newrow[7].replace(' ', '')
        if newrow[7] == '[]':
            base_qual_list_f = []
        else:
            base_qual_list = newrow[7].split(',')
    #       print base_qual_list
            base_qual_list[0] = base_qual_list[0].replace('[', '')
            base_qual_list[-1] = base_qual_list[-1].replace(']', '')
    #       print base_qual_list
            base_qual_list_f = [float(i) for i in base_qual_list]
        return [cell_ftr, read_bases, base_qual_list_f]

    def read_last_line(self, s):
        with open(s, "rb") as f:
            f.seek(-2, 2)
            while (f.read(1) != '\n'):
                f.seek(-2, 1)
            last = f.readline()
        last = last.replace('\n', '')
        pos = int(last.split('\t')[1])
        return pos

    def checkReadPresence(self, single_cell_dict):
        if single_cell_dict.depth == 0:
            return 0
        else:
            return 1

    def Create_Factorial_List(self, max_allele_cnt):
        factorial_list = max_allele_cnt * [None]
        f = math.factorial
        for i in range(max_allele_cnt):
            factorial_list[i] = f(i)
        return factorial_list

    def Create_nCr_mat(self, max_allele_cnt, factorial_list):
        ncr_mat = np.zeros((max_allele_cnt, max_allele_cnt))
        for i in range(max_allele_cnt):
            for j in range(max_allele_cnt):
                ncr_mat[j, i] = factorial_list[j] / \
                    (factorial_list[i] * factorial_list[j - i])
        return ncr_mat

    def CheckAltAllele(self, single_cell_dict):
        """ Check whether the given single cell has alternate allele or not
        """
        alt_c = single_cell_dict.depth - single_cell_dict.refDepth
        if alt_c == 0:
            return 0
        else:
            return 1

    def RefCountString(self, read_base):
        forward_ref_c = read_base.count('.')
        reverse_ref_c = read_base.count(',')
        RefCount = forward_ref_c + reverse_ref_c
        return (forward_ref_c, reverse_ref_c, RefCount)

    def refineBase(self, b):
        b = b.replace(' ', '')
        nb = b.upper()
        return nb

    def copy_list_but_one(self, i_list, index):
        nu_list = []
        len_i_list = len(i_list)
        # print len_i_list
        len_nu_list = len_i_list - 1
        for i in range(index):
            nu_list.append(i_list[i])
        for i in range(index + 1, len_i_list):
            nu_list.append(i_list[i])
        return (nu_list, len_nu_list)

    def update_alt_count(self, alt_allele_count, cell_ftr_dict):
        alt_allele_count[0] += cell_ftr_dict.A_cnt
        alt_allele_count[1] += cell_ftr_dict.T_cnt
        alt_allele_count[2] += cell_ftr_dict.G_cnt
        alt_allele_count[3] += cell_ftr_dict.C_cnt
        return alt_allele_count

    def find_indel(self, string, pattern):
        """ Finds all the occurrences of pattern in string. Removes all 
        occurrences of pattern from the string and returns list of found patterns
        and the new string with patterns removed
        """
        l = [x.group() for x in re.finditer(pattern, string)]
        len_indel = [int(re.split(r'(\d+)', i)[1])
                     for i in l]  # Get the lengths of the patterns
        spans = [i.span() for i in re.finditer(
            pattern, string)]  # Get the spans
        newspan = []  # Change the spans according to the integer length of the pattern
        for i in range(len(len_indel)):
            new_end = spans[i][0] + 1 + len_indel[i] + len(str(len_indel[i]))
            newspan.append((spans[i][0], new_end))
        # Find the instances of the patterns
        final_indel_list = [string[i1:i2] for (i1, i2) in newspan]
        new_string = string
        for i in final_indel_list:
            new_string = new_string.replace(i, '')
        return (final_indel_list, new_string)

    def Alt_count(self, string):
        cp_string = copy.deepcopy(string)
        (ins_list, ins_rmvd_str) = self.find_indel(
            cp_string, '\+[0-9]+[ACGTNacgtn]+')
        ins_count = len(ins_list)
        (del_list, del_ins_rmvd_str) = self. find_indel(
            ins_rmvd_str, '-[0-9]+[ACGTNacgtn]+')
        del_count = len(del_list)
        A_cnt = del_ins_rmvd_str.count('A') + del_ins_rmvd_str.count('a')
        T_cnt = del_ins_rmvd_str.count('T') + del_ins_rmvd_str.count('t')
        G_cnt = del_ins_rmvd_str.count('G') + del_ins_rmvd_str.count('g')
        C_cnt = del_ins_rmvd_str.count('C') + del_ins_rmvd_str.count('c')
        N_cnt = del_ins_rmvd_str.count('N') + del_ins_rmvd_str.count('n')
        return (del_ins_rmvd_str, ins_count, del_count, A_cnt, T_cnt, G_cnt, C_cnt, N_cnt)

    def Count_Start_and_End(self, s):
        end_counts = s.count('$')
        ns = s.replace('$', '')
        start_counts = 0
        i = 0
        fs = ''
        while (i < len(ns)):
            if ns[i] == '^':
                i += 2
                start_counts += 1
            else:
                fs = fs + ns[i]
                i += 1
        return (start_counts, end_counts, fs)

    def Create_base_call_string(self, s, ref):
        """ Removes unwanted characters from s and then replace . and , with ref, 
            finally returns a string that contains all the observed bases
        """
        l = ['.', ',', 'a', 'A', 'c', 'C', 't', 'T', 'g', 'G', '*']
        sn = ''
        for i in s:
            if i in l:
                sn = sn + i
        snn = ''
        for i in sn:
            if i == '.':
                snn = snn + ref
            elif i == ',':
                snn = snn + ref
            elif i == 'a':
                snn = snn + 'A'
            elif i == 'c':
                snn = snn + 'C'
            elif i == 't':
                snn = snn + 'T'
            elif i == 'g':
                snn = snn + 'G'
            elif i == '*':
                snn = snn + ref
            else:
                snn = snn + i
        return snn

    def Get_base_qual_list(self, s):
        """ Returns the base quality scores as a list """
        len_s = len(s)
        base_q_list = len_s * [None]
        base_q_int_list = len_s * [None]
        for i in range(len_s):
            err_p = base_q_tbl.base_q_dict[s[i]]
            err_int = base_q_tbl.base_q_int_dict[s[i]]
            base_q_list[i] = err_p
            base_q_int_list[i] = err_int
        return (base_q_list, base_q_int_list)

    def find_min_list(self, actual_list, flag_list, last_min_loc, last_min_index):
        if ((sum(flag_list) == 1) & (actual_list[last_min_index] == last_min_loc + 1)):
            min_index = last_min_index
            min_loc = actual_list[last_min_index]
        else:
            min_loc = min(actual_list)
            min_index = actual_list.index(min_loc)
        return(min_loc, min_index)

    def feature_row(self, rlist):
        """ Returns the feature row for one cell given the list of input for one cell from pile up file"""
        stat_dict = {}
        if (len(rlist) != 0):
            stat_dict['ref_base'] = self.refineBase(rlist[2])
            stat_dict['depth'] = rlist[3]
            stat_dict['refcount'] = self.RefCountString(rlist[4])
            (del_ins_rmvd_str, ins_count, del_count, A_cnt, T_cnt,
             G_cnt, C_cnt, N_cnt) = self.Alt_count(rlist[4])
            stat_dict['ins_count'] = ins_count
            stat_dict['del_count'] = del_count
            stat_dict['A_count'] = A_cnt
            stat_dict['T_count'] = T_cnt
            stat_dict['G_count'] = G_cnt
            stat_dict['C_count'] = C_cnt
            stat_dict['N_count'] = N_cnt
            stat_dict['base_calls'] = self.Create_base_call_string(
                del_ins_rmvd_str, rlist[2])
            stat_dict['base_qual_list'] = self.Get_base_qual_list(rlist[5])
        else:
            stat_dict['depth'] = 0
            stat_dict['refcount'] = 0
            stat_dict['ref_base'] = ''
            stat_dict['ins_count'] = 0
            stat_dict['del_count'] = 0
            stat_dict['A_count'] = 0
            stat_dict['T_count'] = 0
            stat_dict['G_count'] = 0
            stat_dict['C_count'] = 0
            stat_dict['N_count'] = 0
            stat_dict['base_calls'] = 'NULL'
            stat_dict['base_qual_list'] = []
        return stat_dict

    def calc_strand_bias(self, cell_ftr_pos_list, Alt_count):
        forward_ref_count = 0
        forward_alt_count = 0
        reverse_ref_count = 0
        reverse_alt_count = 0
        for i in range(len(cell_ftr_pos_list)):
            (fr, fa, rr, ra) = cell_ftr_pos_list[i].Store_Strand_Bias_info()
            forward_ref_count += fr
            forward_alt_count += fa
            reverse_ref_count += rr
            reverse_alt_count += ra
        if ((forward_ref_count == 0) & (reverse_ref_count == 0)):
            return (0.0, 1)
        else:
            cont_table = np.array([[forward_ref_count, reverse_ref_count], [
                                  forward_alt_count, reverse_alt_count]])
            (oddsRatio, pval) = stats.fisher_exact(cont_table)
            # if math.isnan(oddsRatio):
            #     for i in range(len(cell_ftr_pos_list)):
            #         print cell_ftr_pos_list[i].start_end_ins_del_rmvd_bases, cell_ftr_pos_list[i].forward_alt_count, cell_ftr_pos_list[i].reverse_alt_count
            #     print cont_table, cell_ftr_pos_list[0].altBase, Alt_count
            return (oddsRatio, pval)

        # if (min(forward_ref_count, reverse_ref_count) == 0):
  #           		refRatio = max(forward_ref_count, reverse_ref_count)
  #           	else:
  #       		refRatio = float(max(forward_ref_count, reverse_ref_count))/min(forward_ref_count, reverse_ref_count)
  #           	if (min(forward_alt_count, reverse_alt_count) == 0):
  #               	altRatio = max(forward_alt_count, reverse_alt_count)
  #           	else:
  #       		altRatio = float(max(forward_alt_count, reverse_alt_count))/min(forward_alt_count, reverse_alt_count)
  #           	if (reverse_ref_count*forward_alt_count == 0):
  #               	R = forward_ref_count*reverse_alt_count
  #           	else:
  #       		R = float(forward_ref_count*reverse_alt_count)/(reverse_ref_count*forward_alt_count)
        # if R == 0:
        # 	oddsRatio = forward_alt_count + forward_ref_count + reverse_ref_count + reverse_alt_count
        # else:
  #       		oddsRatio = R + 1.0/R
  #       	return (refRatio, altRatio, oddsRatio)

    def calc_prior(self, theta, n_cells, flag):
        prior_variant_number = []
        if flag == 1:
            for i in range(0, 2 * n_cells + 1):
                if ((i == 0) | (i == 2 * n_cells)):
                    lst = [1.0 / i for i in range(1, 2 * n_cells)]
                    prob = 0.5 * (1 - (theta * sum(lst)))
                else:
                    prob = theta / i
                prior_variant_number.append(prob)
        elif flag == 2:
            norm_const = 0
            for i in range(1, 2 * n_cells):
                if (i == n_cells):
                    norm_const += 2 * n_cells
                else:
                    norm_const += float(n_cells) / abs(n_cells - i)
            for i in range(1, 2 * n_cells):
                if (i == n_cells):
                    prob = 2 * n_cells * theta / norm_const
                else:
                    prob = ((n_cells * theta) / abs(n_cells - i)) / norm_const
                prior_variant_number.append(prob)
            sp = sum(prior_variant_number)
            p_l0 = 0.5 * (1 - sp)
            p_l2n = 0.5 * (1 - sp)
            prior_variant_number.insert(0, p_l0)
            prior_variant_number.append(p_l2n)
        elif flag == 3:
            for i in range(0, 2 * n_cells + 1):
                prob = 1. / 2 * n_cells + 1
                prior_variant_number.append(prob)
        return prior_variant_number

    def find_max_prob_ratio(self, matrix, dim):
        l_0_prob = matrix.denom_prob_matrix[0, dim[1] - 1]
        if l_0_prob == 0:
            subtracting_max_prob_val = -743.7469
        else:
            subtracting_max_prob_val = math.log(l_0_prob)
        (first_term_max_prob_val, max_prob_allele_count) = max((v, i)
                                                               for i, v in enumerate(matrix.denom_prob_matrix[:, dim[1] - 1]))
        if (first_term_max_prob_val <= 0):
            log_first_term_max_prob_val = -743.7469
        else:
            log_first_term_max_prob_val = math.log(first_term_max_prob_val)
        max_prob_ratio = log_first_term_max_prob_val - subtracting_max_prob_val
        return (max_prob_ratio, max_prob_allele_count)

    def Get_prior_allele_mat(self, read_smpl_count, alt_smpl_count, cell_no_threshold, total_depth, Alt_freq, pe):
        if ((read_smpl_count > cell_no_threshold - 1) & (alt_smpl_count == 1)):
            prior_allele_mat = allele_prior(0.2)
        elif ((read_smpl_count > cell_no_threshold) & (alt_smpl_count == 2) & (total_depth > 30) & (Alt_freq < 0.1)):
            prior_allele_mat = allele_prior(0.1)
        else:
            prior_allele_mat = allele_prior(pe)
        return prior_allele_mat

    def Calc_chr_count(self, barcode):
        AC = 0
        AN = 0
        for c in barcode[1:-1]:
            if c == 'X':
                continue
            else:
                AN += 2
                AC += int(c)
        AF = float(AC) / AN
        return (AC, AF, AN)

    def Calc_Base_Q_Rank_Sum(self, read_supported_cell_list, refBase, altBase):
        ref_read_list = []
        alt_read_list = []
        for cell_ftr_info in read_supported_cell_list:
            for i in range(len(cell_ftr_info.final_bases)):
                if cell_ftr_info.final_bases[i] == refBase:
                    ref_read_list.append(
                        cell_ftr_info.base_qual_int_val_list[i])
                elif cell_ftr_info.final_bases[i] == altBase:
                    alt_read_list.append(
                        cell_ftr_info.base_qual_int_val_list[i])
        return stats.ranksums(alt_read_list, ref_read_list)

    def Calc_qual_depth(self, barcode, All_single_cell_ftrs_list, Qual):
        depth = 0
        for i in range(len(barcode[1:-1])):
            if barcode[i + 1] == 'X':
                continue
            elif barcode[i + 1] == '0':
                continue
            else:
                depth += All_single_cell_ftrs_list[i].depth
        if depth > 0:
            qual_depth = float(Qual) / depth
        else:
            qual_depth = Qual
        return qual_depth

    def Get_BAM_RG(self, bam_file):
        rows = pysam.view("-H", bam_file)
        flag = 0
        for r in rows:
            if r.startswith('@RG'):
                r_l = r.split('\t')
                id = r_l[1].split(':')
                flag = 1
                return id[1]
        if flag == 0:
            bam_id_row = bam_file.split('/')
            bam_id = bam_id_row[-1]
            bam_id = bam_id.replace('..', '')
            bam_id = bam_id.replace('~', '')
            return bam_id

    def Calc_Per_Smpl_Alt_Ref_Ratio(self, total_ref_depth, Alt_count, read_smpl_count, alt_smpl_count):
        if total_ref_depth == 0:
            denom = 1
        else:
            denom = float(total_ref_depth) / read_smpl_count
        num = float(Alt_count) / alt_smpl_count
        return num / denom

    def Consensus_Filter(self, barcode):
        nu_barcode = barcode[1:-1]
        g_count = 0
        for i in nu_barcode:
            if i == '0':
                g_count += 0
            elif i == 'X':
                g_count += 0
            elif i == '1':
                g_count += 1
            elif i == '2':
                g_count += 1
        if g_count > 1:
            return 1
        else:
            return 0
