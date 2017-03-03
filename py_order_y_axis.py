# Standard python template for opening files 
#!/usr/bin/env python

from __future__ import print_function
import argparse
import collections
import csv
import gzip
import magic
import pprint
import re
import operator
import glob
import math

def open_maybe_gzipped(f):
    with open(f) as test_read:
        mimetype = magic.from_buffer(test_read.read(1024), mime=True).decode('utf-8')
    if mimetype == 'application/x-gzip':
        f = gzip.open(f, mode='rt') # Python3: , encoding='utf-8')
    else:
        f = open(f, 'rt')
    return f

if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('datfile', help="""The data file to format""")
    parser.add_argument('orderfile', help="""The file with cell types ranked according to the order.""")
    args = parser.parse_args()
    # print("trait\tenhancer_type\tcell\tthreshold\toverlap\texpected_overlap\tpval\tcorrected_pval\tlog2_enrichment\tsigned_minuslog_pval\tsigned_minuslog_corrected_pval\tfilt_log2_enrich\tfilt_signed_minuslog_corrected_pval")

    with open_maybe_gzipped(args.datfile) as d_file, open_maybe_gzipped(args.orderfile) as o_file :
        dat_file = csv.reader(d_file, delimiter='\t', dialect='excel-tab')
        order_file = csv.reader(o_file, delimiter='\t', dialect='excel-tab')
        for stuff in dat_file:
            o_file.seek(0)
            for line in order_file:
                # print(line)
                # print(stuff)
                if str(stuff[1]) == str(line[0]):
                    stuff.append(str(line[1]))
                elif str(stuff[1]) == 'cell':
                    stuff.append("y_order")
                    break

                
            print(*stuff, sep='\t')
