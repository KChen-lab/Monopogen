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

from alleles_prior import allele_prior
from utils import Utils_Functions
import math
import copy

U = Utils_Functions()


class Single_Cell_Ftrs_Pos:

    # Constructor takes the list of info for the current position as input
    # ([contig, loc, ref, depth, primary_bases, base_q])
    def __init__(self, refBase, current_pos_info_list):
        if (int(current_pos_info_list[0]) == 0):
            self.depth = 0
            self.refDepth = 0
        else:
            self.refBase = refBase
            self.depth = int(current_pos_info_list[0])
            self.primary_bases = current_pos_info_list[1]
            self.base_q = current_pos_info_list[2]
            (self.forward_ref_count, self.reverse_ref_count,
             self.refDepth) = U.RefCountString(self.primary_bases)

    # Remove the insertions and deletions from the primary_bases and also
    # count the number of insertions and deletions
    def Get_Ins_Del_rmvd_bases(self):
        if ((self.primary_bases.count('+') + self.primary_bases.count('-')) == 0):
            self.ins_count = 0
            self.del_count = 0
            self.ins_list = []
            self.del_list = []
            self.ins_del_rmvd_bases = self.primary_bases
        else:
            cp_primary_bases = copy.deepcopy(self.primary_bases)
            (self.ins_list, ins_rmvd_bases) = U.find_indel(
                cp_primary_bases, '\+[0-9]+[ACGTNacgtn]+')
            self.ins_count = len(self.ins_list)
            (self.del_list, self.ins_del_rmvd_bases) = U.find_indel(
                ins_rmvd_bases, '-[0-9]+[ACGTNacgtn]+')
            self.del_count = len(self.del_list)
        return 0

    # Function that calculates the base calling errors from base qual scores
    def Get_Base_Qual_Vals(self):
        (self.base_qual_val_list,
         self.base_qual_int_val_list) = U.Get_base_qual_list(self.base_q)
        return 0

    # After removal of insertions and deletions we create the base call string
    # to be used by the model
    def Get_Base_Calls(self, ref):
        (self.start_read_counts, self.end_read_counts,
         self.start_end_ins_del_rmvd_bases) = U.Count_Start_and_End(self.ins_del_rmvd_bases)
        self.final_bases = U.Create_base_call_string(
            self.start_end_ins_del_rmvd_bases, ref)
        return 0

    def Get_base_call_string_nd_quals(self, ref):
        self.Get_Ins_Del_rmvd_bases()
        self.Get_Base_Qual_Vals()
        self.Get_Base_Calls(ref)
        # print self.final_bases, self.base_qual_val_list
        return 0

    # Function that calculates the numbers of alternate alleles in the
    # ins_del_rmvd_bases
    def Get_Alt_Allele_Count(self):
        self.A_cnt = self.start_end_ins_del_rmvd_bases.count(
            'A') + self.start_end_ins_del_rmvd_bases.count('a')
        self.C_cnt = self.start_end_ins_del_rmvd_bases.count(
            'C') + self.start_end_ins_del_rmvd_bases.count('c')
        self.G_cnt = self.start_end_ins_del_rmvd_bases.count(
            'G') + self.start_end_ins_del_rmvd_bases.count('g')
        self.T_cnt = self.start_end_ins_del_rmvd_bases.count(
            'T') + self.start_end_ins_del_rmvd_bases.count('t')
        return 0

    # Saves the cell index
    def Set_Cell_Index(self, index):
        self.cell_index = index
        return 0

    # Store altBase, Alt_freq, prior_allele_mat
    def Store_Addl_Info(self, refBase, altBase, Alt_freq, prior_allele_mat):
        self.altBase = altBase
        self.Alt_freq = Alt_freq
        self.prior_allele_mat = prior_allele_mat
        self.refBase = refBase
        # self.altBase          = U.refineBase(self.altBase)
        return 0

    def Store_Strand_Bias_info(self):
        # self.forward_ref_count = self.start_end_ins_del_rmvd_bases.count('.')
        # self.reverse_ref_count = self.start_end_ins_del_rmvd_bases.count(',')
        self.lowercase_alt_base = self.altBase.lower()
        self.forward_alt_count = self.start_end_ins_del_rmvd_bases.count(
            self.altBase)
        self.reverse_alt_count = self.start_end_ins_del_rmvd_bases.count(
            self.lowercase_alt_base)
        self.altcount = self.forward_alt_count + self.reverse_alt_count
        return (self.forward_ref_count, self.forward_alt_count, self.reverse_ref_count, self.reverse_alt_count)

    # Refing g
    def refineG(self, g):
        if g == 'CA':
            g = 'AC'
        elif g == 'GA':
            g = 'AG'
        elif g == 'TA':
            g = 'AT'
        elif g == 'GC':
            g = 'CG'
        elif g == 'TC':
            g = 'CT'
        elif g == 'TG':
            g = 'GT'
        return g

    # Calculate probability of data given genotype gt
    def Calc_Prob_gt(self, gt, max_depth):
        val = 1.0
        ub = min(len(self.base_qual_val_list),
                 len(self.final_bases), max_depth)
        for i in range(ub):
            curr_base = self.final_bases[i]
            # curr_base     = U.refineBase(curr_base)
            curr_base_key = (gt, curr_base)
            curr_err = self.base_qual_val_list[i]
            prob_i = self.prior_allele_mat.getValue(curr_base_key)
            prob = curr_err * (1 - prob_i) / 3 + (1 - curr_err) * prob_i
            val = val * prob
        return val

    # Function to calculate likelihood of homo ref genotype for ub = 50
    def Prob_Reads_Given_Genotype_homo_ref_50d(self, g, ub):
        key_curr_base = (g, self.refBase)
        curr_base_genotype_prob = self.prior_allele_mat.getValue(key_curr_base)
        complement_curr_base_genotype_prob = 1 - curr_base_genotype_prob
        probability = 1.0
        for i in range(ub):
            curr_base = self.final_bases[i]
            curr_err = self.base_qual_val_list[i]
            if (curr_base == self.refBase):
                prob_i = curr_err * complement_curr_base_genotype_prob / \
                    3 + (1 - curr_err) * curr_base_genotype_prob
            else:
                prob_i = curr_err * curr_base_genotype_prob + \
                    (1 - curr_err) * complement_curr_base_genotype_prob / 3
            probability = probability * prob_i
        self.cell_prob_0 = probability
        self.cell_prob_0_50d = probability
        return self.cell_prob_0

    # Function to calculate likelihood of homo ref genotype for ub = all
    def Prob_Reads_Given_Genotype_homo_ref_all(self, g, ub):
        key_curr_base = (g, self.refBase)
        curr_base_genotype_prob = self.prior_allele_mat.getValue(key_curr_base)
        complement_curr_base_genotype_prob = 1 - curr_base_genotype_prob
        probability = 1.0
        for i in range(100):
            curr_base = self.final_bases[i]
            curr_err = self.base_qual_val_list[i]
            if (curr_base == self.refBase):
                prob_i = curr_err * complement_curr_base_genotype_prob / \
                    3 + (1 - curr_err) * curr_base_genotype_prob
            else:
                prob_i = curr_err * curr_base_genotype_prob + \
                    (1 - curr_err) * complement_curr_base_genotype_prob / 3
            probability = probability * prob_i
        self.cell_prob_0_50d = probability
        for i in range(100, ub):
            curr_base = self.final_bases[i]
            curr_err = self.base_qual_val_list[i]
            if (curr_base == self.refBase):
                prob_i = curr_err * complement_curr_base_genotype_prob / \
                    3 + (1 - curr_err) * curr_base_genotype_prob
            else:
                prob_i = curr_err * curr_base_genotype_prob + \
                    (1 - curr_err) * complement_curr_base_genotype_prob / 3
            probability = probability * prob_i
        self.cell_prob_0 = probability

        return self.cell_prob_0

    # Function to calculate likelihood of homo nonref genotype for ub = 50
    def Prob_Reads_Given_Genotype_homo_nonref_50d(self, g, ub):
        key_curr_base = (g, self.altBase)
        curr_base_genotype_prob = self.prior_allele_mat.getValue(key_curr_base)
        complement_curr_base_genotype_prob = 1 - curr_base_genotype_prob
        probability = 1.0
        for i in range(ub):
            curr_base = self.final_bases[i]
            curr_err = self.base_qual_val_list[i]
            if (curr_base == self.altBase):
                prob_i = curr_err * complement_curr_base_genotype_prob / \
                    3 + (1 - curr_err) * curr_base_genotype_prob
            else:
                prob_i = curr_err * curr_base_genotype_prob + \
                    (1 - curr_err) * complement_curr_base_genotype_prob / 3
            probability = probability * prob_i
        self.cell_prob_2 = probability
        self.cell_prob_2_50d = probability
        return self.cell_prob_2

    # Function to calculate likelihood of homo nonref genotype for ub = all
    def Prob_Reads_Given_Genotype_homo_nonref_all(self, g, ub):
        key_curr_base = (g, self.altBase)
        curr_base_genotype_prob = self.prior_allele_mat.getValue(key_curr_base)
        complement_curr_base_genotype_prob = 1 - curr_base_genotype_prob
        probability = 1.0
        for i in range(100):
            curr_base = self.final_bases[i]
            curr_err = self.base_qual_val_list[i]
            if (curr_base == self.altBase):
                prob_i = curr_err * complement_curr_base_genotype_prob / \
                    3 + (1 - curr_err) * curr_base_genotype_prob
            else:
                prob_i = curr_err * curr_base_genotype_prob + \
                    (1 - curr_err) * complement_curr_base_genotype_prob / 3
            probability = probability * prob_i
        self.cell_prob_2_50d = probability
        for i in range(100, ub):
            curr_base = self.final_bases[i]
            curr_err = self.base_qual_val_list[i]
            if (curr_base == self.altBase):
                prob_i = curr_err * complement_curr_base_genotype_prob / \
                    3 + (1 - curr_err) * curr_base_genotype_prob
            else:
                prob_i = curr_err * curr_base_genotype_prob + \
                    (1 - curr_err) * complement_curr_base_genotype_prob / 3
            probability = probability * prob_i
        self.cell_prob_2 = probability

        return self.cell_prob_2

    # Helper function to calculate beyond 50 bases
    def Calc_Prob_gt_beyond_50(self, g, ub, prob_1_50d):
        val = prob_1_50d
        for i in range(100, ub):
            curr_base = self.final_bases[i]
            curr_base_key = (g, curr_base)
            curr_err = self.base_qual_val_list[i]
            prob_i = self.prior_allele_mat.getValue(curr_base_key)
            prob = curr_err * (1 - prob_i) / 3 + (1 - curr_err) * prob_i
            val = val * prob
        return val

    # Function to calculate likelihood of hetero genotype for ub = 50
    def Prob_Reads_Given_Genotype_hetero_50d(self, g, ub, pad):
        prob_0 = self.cell_prob_0_50d
        prob_2 = self.cell_prob_2_50d
        prob_1 = self.Calc_Prob_gt(g, ub)
        pad_c = pad / 2
        probability = (1 - pad) * prob_1 + pad_c * prob_0 + pad_c * prob_2
        self.cell_prob_1 = probability
        self.cell_prob_1_50d = probability
        return self.cell_prob_1

    # Function to calculate likelihood of hetero genotype for ub = all
    def Prob_Reads_Given_Genotype_hetero_all(self, g, ub, pad):
        prob_0_50d = self.cell_prob_0_50d
        prob_2_50d = self.cell_prob_2_50d
        prob_1_50d = self.Calc_Prob_gt(g, 100)
        # prob_1_50d = self.Calc_Prob_gt(g, 50)
        pad_c = pad / 2
        probability = (1 - pad) * prob_1_50d + pad_c * \
            prob_0_50d + pad_c * prob_2_50d
        self.cell_prob_1_50d = probability
        prob_0 = self.cell_prob_0
        prob_2 = self.cell_prob_2
        prob_1 = self.Calc_Prob_gt_beyond_50(g, ub, prob_1_50d)
        all_prob = (1 - pad) * prob_1 + pad_c * prob_0 + pad_c * prob_2
        self.cell_prob_1 = all_prob
        return self.cell_prob_1

    def Prob_Reads_Given_Genotype(self, genotype_flag, max_depth, pad):
        # print "called cell and gt", self.cell_index, genotype_flag
        if (self.altBase == ''):
            if (genotype_flag != 0):
                self.cell_prob_2 = 0
                self.cell_prob_1 = 0

                return 0.0
            else:
                g = self.refBase + self.refBase
        else:
            if (genotype_flag == 0):
                g = self.refBase + self.refBase
            elif (genotype_flag == 2):
                g = self.altBase + self.altBase
            else:
                g = self.refBase + self.altBase
                g = self.refineG(g)
        # print g
        # print "called cell and gt", self.cell_index, genotype_flag,
        # self.final_bases, g
        ub = min(len(self.base_qual_val_list),
                 len(self.final_bases), max_depth)
        if (ub <= 100):
            if (genotype_flag == 0):
                probability = self.Prob_Reads_Given_Genotype_homo_ref_50d(
                    g, ub)
                # print ub, genotype_flag, self.cell_prob_0
                return probability
            elif (genotype_flag == 2):
                probability = self.Prob_Reads_Given_Genotype_homo_nonref_50d(
                    g, ub)
                # print ub, genotype_flag, self.cell_prob_2
                return probability
            elif (genotype_flag == 1):
                probability = self.Prob_Reads_Given_Genotype_hetero_50d(
                    g, ub, pad)
                # print ub, genotype_flag, self.cell_prob_1
                return probability
        else:
            if (genotype_flag == 0):
                probability = self.Prob_Reads_Given_Genotype_homo_ref_all(
                    g, ub)
                # print ub, genotype_flag, probability
                return probability
            elif (genotype_flag == 2):
                probability = self.Prob_Reads_Given_Genotype_homo_nonref_all(
                    g, ub)
                # print ub, genotype_flag, probability
                return probability
            elif (genotype_flag == 1):
                probability = self.Prob_Reads_Given_Genotype_hetero_all(
                    g, ub, pad)
                # print ub, genotype_flag, probability
                return probability

    def Prob_Reads_Given_Genotype_50d(self, genotype_flag):
        if (genotype_flag == 0):
            self.cell_prob_0 = self.cell_prob_0_50d
            return self.cell_prob_0_50d
        elif (genotype_flag == 2):
            self.cell_prob_2 = self.cell_prob_2_50d
            return self.cell_prob_2_50d
        elif (genotype_flag == 1):
            self.cell_prob_1 = self.cell_prob_1_50d
            return self.cell_prob_1_50d

    def Prob_Reads_Given_Genotype_Genotyping(self, gt_flag):
        if gt_flag == 0:
            return self.cell_prob_0
        elif gt_flag == 1:
            return self.cell_prob_1
        else:
            return self.cell_prob_2
