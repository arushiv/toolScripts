# Input a matrix with gene name and ESI scores for tissues to partition by fractions specified 
#!/usr/bin/env python

import argparse
import collections
import csv
import gzip
import magic
import pprint
import re
import operator
import glob
import pandas
import numpy as np

def open_maybe_gzipped(f):
    with open(f) as test_read:
        mimetype = magic.from_buffer(test_read.read(1024), mime=True).decode('utf-8')
    if mimetype == 'application/x-gzip':
        f = gzip.open(f, mode='rt') # Python3: , encoding='utf-8')
    else:
        f = open(f, 'rt')
    return f

class bestEqtl(object):
    def __init__(self):
        pass
    
    def subset_by_pval_threshold(self, d, pval, gene, threshold):
        self.pval = pval
        self.gene = gene
        dsubset = d[d[pval] <= threshold]
        return dsubset
        
    def aggregate_by_best_pval(self, dsubset):
        # dnew = dsubset.groupby(self.header).min()[self.pval]
        # dnew = dsubset.groupby(self.gene).min()[self.pval]
        dnew = dsubset.loc[dsubset.groupby(self.gene)[self.pval].idxmin()]
        return dnew
        # dnew.columns = dnew.columns.droplevel(0)
        # print pandas.merge(dsubset, dnew, how="right", on = self.gene)
        # print pandas.merge(dnew, dsubset)
    # return d.sort_values(by = name).dropna()
    # dsort = d[d[name] != 'NA'].sort_values(by = name)
    # return dsort

    def print_to_file(self, dnew, outfile):
        # self.header.append(self.pval)
        dnew.to_csv(outfile, sep='\t', index=False)

if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('datfile', type=str, help="""The tab separated data file with at least gene name, and p values for eQTLs.""")
    parser.add_argument('-pn','--pval_column', type=str, nargs='?', default='pvalue', help="""Column name for pvalue.""")
    parser.add_argument('-gn','--gene_column', type=str, nargs='?', default='gene', help="""Column name for gene identifiers.""")
    parser.add_argument('-t','--threshold', type=float, nargs='?', default=0.05, help="""p value threshold. Default=0.05""")
    parser.add_argument('outfile', type=str, help="""output file name.""")
    
    args = parser.parse_args()
    
    d = pandas.read_csv(args.datfile, sep='\t')
    gene = args.gene_column
    pval = args.pval_column
    threshold = args.threshold
    outfile = args.outfile

    f = bestEqtl()
    f.print_to_file(f.aggregate_by_best_pval(f.subset_by_pval_threshold(d, pval, gene, threshold)), outfile)
