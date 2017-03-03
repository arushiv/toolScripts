# Match gene name on provided column, separate into files as per quantiles.
# Ex. binQuintile file at /lab/arushiv/chromatin/2015_12_08_islet_eQTL/eqtl_binnedByIsletGenes/bins.protein_coding.binQuintile.1.byGenes contains gene names that are in Quintile 1
#!/usr/bin/env python


#from __future__ import print_function
import argparse
import collections
import csv
import gzip
import magic
import pprint
import re
import operator
import glob

def open_maybe_gzipped(f):
    with open(f) as test_read:
        mimetype = magic.from_buffer(test_read.read(1024), mime=True).decode('utf-8')
    if mimetype == 'application/x-gzip':
        f = gzip.open(f, mode='rt') # Python3: , encoding='utf-8')
    else:
        f = open(f, 'rt')
    return f

def splitByQuantile(quantilef, eqtlf, column_gene, qf, outputfile):
    with open_maybe_gzipped(eqtlf) as efile:
        eqtlfile = csv.reader(efile, delimiter='\t', dialect='excel-tab')

        with open_maybe_gzipped(quantilef) as qfile:
            quantilefile = csv.reader(qfile, delimiter='\t', dialect='excel-tab')

            with open(outputfile, 'w') as fout:
                for line in eqtlfile:
                    for stuff in quantilefile:
                        if line[column_gene - 1] == stuff[qf - 1]:   ### If gene name format in both files is same
                        # if (line[column_gene - 1]).split('.')[0] == stuff[qf - 1]:      ### If gene name format in quantile file is not of the format ENSG00XYZW.K
                            fout.write('\t'.join(line))
                            fout.write('\n')
                    qfile.seek(0)
                    
if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('quantilefile', help="""The file with genes for a certain quantile.""")
    parser.add_argument('eqtlfile', help="""The file with data to be searched and split.""")
    parser.add_argument('-f', '--column_gene', type=int, default=1, help="""The column number in the eqtl file which contains the gene name.""")
    parser.add_argument('-qf', '--column_quantileFile', type=int, default=1, help="""The column number in quantile file for the gene name.""")
    parser.add_argument('outputfile', help="""Name of the output file.""")

    args = parser.parse_args()

    quantilef = args.quantilefile
    eqtlf = args.eqtlfile
    column_gene = args.column_gene
    outputfile = args.outputfile
    qf = args.column_quantileFile
    
    splitByQuantile(quantilef, eqtlf, column_gene, qf, outputfile)
