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

def checkDuplicateGenes(d, gene):
    # print d.iloc[:,0]
    # print d.duplicated(subset = 'gene')
    if d.set_index(gene).index.get_duplicates():
        print "Duplicate gene names exist! Check input file"

def split_dataframe(dsort, quantiles):
    return np.array_split(dsort, quantiles)

def savefiles(listdf, gene, name, outname):
    for i in range(len(listdf)):
        namefile = outname + "." + str(i + 1) + ".txt"
        header = [gene, name]
        listdf[i].to_csv(namefile, sep='\t', index=False, columns=header)
    
def sort_by_esi(d, name):
    return d.sort_values(by = name).dropna()
    # dsort = d[d[name] != 'NA'].sort_values(by = name)
    # return dsort


if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('esi_matrix', type=str, help="""The tab separated matrix with genes (column 1) and ESI score in one tissue across subsequent columns.""")
    parser.add_argument('-n','--tissue_name', type=str, help="""Column name for tissue to sort ESI scores for""")
    parser.add_argument('-gn','--gene_column', type=str, nargs='?', default='_gene', help="""Column name for gene identifiers""")
    parser.add_argument('-f','--quantiles', type=int, nargs='?', default=5, help="""Number of quantiles to partition genes in""")
    parser.add_argument('outFileIdentifier', type=str, help="""output files will be named as outFileIdentifier.quantile.txt""")
    
    args = parser.parse_args()
    
    d = pandas.read_csv(args.esi_matrix, sep='\t')
    name = args.tissue_name
    quantiles = args.quantiles
    outname = args.outFileIdentifier
    gene = args.gene_column

    checkDuplicateGenes(d, gene)
    savefiles(split_dataframe(sort_by_esi(d, name), quantiles), gene, name, outname) 
