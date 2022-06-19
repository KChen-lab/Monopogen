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

import argparse
import os
import sys
import gzip
from subprocess import Popen, PIPE, check_call
import re
import time
from glob import glob


VCF_meta_template = """##fileformat=VCFv4.1
##fileDate={_t.tm_year}-{_t.tm_mon}-{_t.tm_mday}
##source=MonoVar
{_d.FILTER_META}{_d.INFO_META}{_d.FORMAT_META}{_d.REF_META}#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{_d.FILES_META}
"""

VCF_record_template = "{_r.CHROM}\t{_r.POS}\t{_r.ID}\t{_r.REF}\t{_r.ALT}\t{_r.QUAL}\t{_r.FILTER}\t{_r.INFO}\t{_r.FORMAT}\t{_r.PASSCODE}\n"


class VCFRecord():

    def __init__(self):

        self.info = {}


class VCFDocument():

    def __init__(self, outf):

        self.time = time.ctime()

        self.info_fields = []
        self.filter_fields = []
        self.format_fields = []
        self.outf = outf

    def populate_fields(self, bam_id_list):
        self.filter_fields.append(('LowQual', 'Low quality'))
        self.format_fields.append(
            ('AD', '.', 'Integer', 'Allelic depths for the ref and alt alleles in the order listed'))
        self.format_fields.append(
            ('DP', '1', 'Integer', 'Approximate read depth (reads with MQ=255 or with bad mates are filtered)'))
        self.format_fields.append(('GQ', '1', 'Integer', 'Genotype Quality'))
        self.format_fields.append(('GT', '1', 'String', 'Genotype'))
        self.format_fields.append(
            ('PL', 'G', 'Integer', 'Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification'))
        self.info_fields.append(
            ('AC', 'A', 'Integer', 'Allele count in genotypes, for each ALT allele, in the same order as listed'))
        self.info_fields.append(
            ('AF', 'A', 'Float', 'Allele Frequency, for each ALT allele, in the same order as listed'))
        self.info_fields.append(
            ('AN', '1', 'Integer', 'Total number of alleles in called genotypes'))
        self.info_fields.append(
            ('BaseQRankSum', '1', 'Float', 'Z-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities'))
        self.info_fields.append(
            ('DP', '1', 'Integer', 'Approximate read depth; some reads may have been filtered'))
        self.info_fields.append(
            ('QD', '1', 'Float', 'Variant Confidence/Quality by Depth'))
        self.info_fields.append(
            ('SOR', '1', 'Float', 'Symmetric Odds Ratio of 2x2 contingency table to detect strand bias'))
        self.info_fields.append(
            ('MPR', '1', 'Float', 'Log Odds Ratio of maximum value of probability of observing non-ref allele to the probability of observing zero non-ref allele'))
        self.info_fields.append(
            ('PSARR', '1', 'Float', 'Ratio of per-sample Alt allele supporting reads to Ref allele supporting reads'))
        self.files_list = bam_id_list

    def populate_reference(self, ref_file):
        self.ref_file = ref_file

    def print_header(self):

        self.FILTER_META = ''.join(
            '##FILTER=<ID=%s,Description="%s">\n' % _ for _ in self.filter_fields)
        self.FORMAT_META = ''.join(
            '##FORMAT=<ID=%s,Number=%s,Type=%s,Description="%s">\n' % _ for _ in self.format_fields)
        self.INFO_META = ''.join(
            '##INFO=<ID=%s,Number=%s,Type=%s,Description="%s">\n' % _ for _ in self.info_fields)
        self.FILES_META = '\t'.join(self.files_list)
        self.REF_META = '##reference=file:{0}\n'.format(self.ref_file)
        self.outf.write(VCF_meta_template.format(_d=self, _t=time.localtime()))

    def print_record(self, record):

        record.INFO = ';'.join("%s=%s" % (
            _[0], str(record.info[_[0]])) for _ in self.info_fields if _[0] in record.info)
        self.outf.write(VCF_record_template.format(_r=record))

    def print_my_record(self, record):
        self.outf.write(VCF_record_template.format(_r=record))

    def close(self):

        self.outf.close()


class VRecord:

    def __init__(self, chrm, pos):
        self.CHROM = chrm
        self.POS = pos

    def get6fields(self, ref, alt, id, qual, filter_i, info):
        self.ID = id
        self.REF = ref
        self.ALT = alt
        self.QUAL = str(qual)
        self.FILTER = filter_i
        self.INFO = info

    def get_passcode(self, barcode):
        self.PASSCODE = barcode

    # def __init__(self, chrm, pos):
    #     self.chrm = chrm
    #     self.pos = int(pos)

    def format_vcf(self, samples):
        self.FORMAT = "\t".join(samples)
        # return "%s\t%d\t%s\t%s\n" % (self.chrm, self.pos, "\t".join(self.text2to8), "\t".join([self.data[_] for _ in samples]))
        # return "%s\t%d\t%s\t%s\n" % (self.chrm, self.pos,
        # "\t".join(text2to8), "\t".join(samples))


class VCF:

    def __init__(self, fn):
        self.fn = fn
        self.header = ''
        self.chrmline = ''

    def open(self):
        if self.fn.endswith(".gz"):
            self.fh = gzip.open(self.fn, "r")
        else:
            self.fh = open(self.fn, "r")

    def close(self):
        self.fh.close()

    def read_header(self):

        while (True):
            line = self.fh.readline()
            if not line:
                break

            if line.startswith("#CHROM"):
                self.samplelist = line.strip().split()[9:]
                self.chrmline = '\t'.join(line.split()[:9])
                break
            else:
                self.header += line

    def format_header(self):
        return self.header + self.chrmline + '\t' + '\t'.join(self.samplelist)

    def format1(self, r):

        s = r.chrm
        s += '\t%d\t' % r.pos
        s += '\t'.join(r.text2to8)
        s += '\t'
        s += '\t'.join([r.data[_] for _ in self.samplelist])

        return s

    def read1(self):

        while (True):
            line = self.fh.readline()
            if not line:
                break
            if line[0] == '#':
                continue
            pair = line.split("\t")
            r = VRecord(pair[0], pair[1])
            r.text2to8 = pair[2:9]
            r.id = pair[2]
            r.data = dict(
                zip(self.samplelist, [_.split(":")[0] for _ in pair[9:]]))
            yield r

    def fetch_region(self, chrm, beg, end):

        for line in Popen(["tabix", self.fn, "%s:%d-%d" % (chrm, beg, end)], stdout=PIPE).stdout:
            pair = line.split("\t")
            r = VRecord(pair[0], pair[1])
            r.text2to8 = pair[2:9]
            r.data = dict(
                zip(self.samplelist, [_.split(":")[0] for _ in pair[9:]]))
            yield r


def main_compare(args):

    vcf1 = VCF(sys.argv[1])
    vcf2 = VCF(sys.argv[2])

    vcf1.open()
    vcf2.open()

    vcf1.read_header()
    vcf2.read_header()

    intersample = set(vcf2.samplelist) & set(vcf1.samplelist)
    print (len(intersample) +  " overlapping samples")

    cs2g2 = {}
    for r in vcf2.read1():
        for s in intersample:
            cs2g2[(r.id, s)] = r.data[s]

    cs2g1 = {}
    for r in vcf1.read1():
        for s in intersample:
            cs2g1[(r.id, s)] = r.data[s]

    print (len(cs2g1) +  " call-sample pairs in vcf1")
    print (len(cs2g2) +  " call-sample pairs in vcf2")

    overlap = set(cs2g1.keys()) & set(cs2g2.keys())
    print (len(overlap) + "overlapping call-sample pairs between the two vcfs")

    vars = ["./.", "0/0", "0/1", "1/1"]
    #print 'vcf1\\vcf2', '\t'.join(vars)
    for var1 in vars:
        sys.stdout.write(var1)
        for var2 in vars:
            sys.stdout.write('\t%d' % len(
                [_ for _ in overlap if cs2g1[_] == var1 and cs2g2[_] == var2]))
        sys.stdout.write('\n')


def main_filter(args):

    vcf = VCF(args.v)
    vcf.open()
    vcf.read_header()

    cids = set()
    with open(args.cid) as fh:
        for line in fh:
            pair = line.split('\t')
            cids.add(pair[args.cidcol - 1])

    print (vcf.format_header())
    for r in vcf.read1():
        if r.id in cids:
            print (vcf.format1(r))

    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='vcf tool')
    subparsers = parser.add_subparsers()

    psr_compare = subparsers.add_parser("compare", help=""" compare vcfs """)
    psr_compare.add_argument('-v1', help='VCF file 1')
    psr_compare.add_argument('-v2', help='VCF file 2')
    psr_compare.set_defaults(func=main_compare)

    psr_filter = subparsers.add_parser("filter", help=""" filter vcf """)
    psr_filter.add_argument('-v', help='VCF file')
    psr_filter.add_argument('--cid', help='call id list')
    psr_filter.add_argument('--cidcol', type=int, default=1,
                            help='call id column index (1-based) [1]')
    psr_filter.set_defaults(func=main_filter)

    args = parser.parse_args()
    args.func(args)
