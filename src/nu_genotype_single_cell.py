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
from Single_Cell_Ftrs_Pos import Single_Cell_Ftrs_Pos
from genotype_prob_mat import Genotype_Prob_matrix
import sys
import math


class Single_cell_genotype_records:

    def __init__(self, current_cell_ftr_list, other_cells_ftr_list, curr_cell_residual_mat, other_cell_list_len, n_cells, prior_allele_mat, prior_variant_allele):
        self.current_cell_ftr_list = current_cell_ftr_list
        self.other_cells_ftr_list = other_cells_ftr_list
        self.residual_mat = curr_cell_residual_mat
        self.other_cell_list_len = other_cell_list_len
        # self.j			           = gt_flag
        self.n_cells = n_cells
        self.prior_allele_mat = prior_allele_mat
        self.prior_variant_allele = prior_variant_allele
        self.p1 = []
        self.p2 = []
    # def current_cell_genotype_prob(self):
    # return
    # self.current_cell_ftr_list.Prob_Reads_Given_Genotype_Genotyping(self.j)

    def current_cell_genotype_prob(self, gt_flag):
        return self.current_cell_ftr_list.Prob_Reads_Given_Genotype_Genotyping(gt_flag)

    def nCr(self, n, r, factorial_list):
        comb_denom = n - r
        number = factorial_list[n] / \
            (factorial_list[r] * factorial_list[comb_denom])
        return number

    # def find_coeff(self, n_cells, l, j, factorial_list):
    # 	if (j > l):
    # 		return 0
    # 	else:
    # 		num = self.nCr(l, j, factorial_list) * self.nCr(2*n_cells - l, 2 - j, factorial_list)
    # 		denom = self.nCr(2*n_cells, 2, factorial_list)
    # 		return float(num)/denom

    def find_coeff(self, n_cells, l, j, nCr_matrix):
        if (j > l):
            return 0
        else:
            num = nCr_matrix[l, j] * nCr_matrix[2 * n_cells - l, 2 - j]
            denom = nCr_matrix[2 * n_cells, 2]
            return float(num) / denom

    def other_cells_genotype_prob(self, gt_flag, nCr_matrix):
        if self.other_cell_list_len == 0:
            return 1
        else:
            prob = 0.0
            for l in range(0, self.residual_mat.dim[0]):
                coeff = self.find_coeff(self.n_cells, l, gt_flag, nCr_matrix)
                prob = prob + coeff * self.residual_mat.denom_prob_matrix[
                    l, self.residual_mat.dim[1] - 1] * self.prior_variant_allele[l + gt_flag]
            return prob

    # def other_cells_genotype_prob(self, max_depth):
    # 	if len(self.other_cells_ftr_list) == 0:
    # 		return 1
    # 	else:
    # 		prob = 0.0
    # 		for l in range(0, self.residual_mat.dim[0]):
    # 			coeff = self.find_coeff(self.n_cells, l, self.j)
    # 			prob = prob + coeff * self.residual_mat.denom_prob_matrix[l, self.residual_mat.dim[1]-1] * self.prior_variant_allele[l+self.j]
    # 		return prob

    def find_genotype_prob(self, gt_flag, nCr_matrix):
        p1 = self.current_cell_genotype_prob(gt_flag)
        p2 = self.other_cells_genotype_prob(gt_flag, nCr_matrix)
        # print p1, p2, gt_flag
#		self.p1.append(p1)
#		self.p2.append(p2)
        if p2 == 0:
            p2 = 1e-322
        if p1 == 0:
            p1 = 1e-322
        return p1 * p2
