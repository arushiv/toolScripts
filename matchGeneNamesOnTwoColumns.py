# Match gene name on provided column, add matching entries into separate files.
# Ex. binQuintile file at /lab/arushiv/chromatin/2015_12_08_islet_eQTL/eqtl_binnedByIsletGenes/bins.protein_coding.binQuintile.1.byGenes contains gene names that are in Quintile 1
#!/usr/bin/env python


#from __future__ import print_function
import argparse
import collections
import csv
import gzip
import pprint
import re
import operator
import glob
import pandas


def splitByQuantile_full(file1, file2, column_1, column_2):
    colname1 = "%s_x"%column_1
    colname2 = "%s_y"%column_2
    df = pandas.merge(file1, file2, how='inner', left_on = column_1, right_on = column_2).drop([colname1, colname2],1)
    return df

def splitByQuantile_partial(file1, file2, column_1, column_2):
    file1['keyColumn'] = file1.iloc[:, column_1].str.split('.').apply(lambda x: x[0])

    file2['keyColumn'] = file2.iloc[:, column_2].str.split('.').apply(lambda x: x[0])
    
    df = pandas.merge(file1, file2, how='inner', on='keyColumn').drop('keyColumn', 1)
    return df

def writeToFile(df, outputfile):
    df.to_csv(outputfile, sep='\t', index=False)

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser = argparse.ArgumentParser(description='Perform inner join on a specified column in two files')
    parser.add_argument('file1', help="""The first file with column to match on, tab delimited, no header""")
    parser.add_argument('file2', help="""The second file with column to match on, tab delimited, no header""")
    parser.add_argument('-c1', '--column_1', type=int, default=1, help="""The column number in the first file.""")
    parser.add_argument('-c2', '--column_2', type=int, default=1, help="""The column number in the second file.""")
    parser.add_argument('outputfile', help="""Name of the output file.""")
    parser.add_argument('-matchType', type=str, default="full", help="""'full' if gene name in two files match exactly. 'partial' if names in both files do not contain .x extensions. Default = 'full'""")

    args = parser.parse_args()

    file1 = pandas.read_csv(args.file1, sep='\t', header = None, dtype=str)
    file2 = pandas.read_csv(args.file2, sep='\t', header = None, dtype=str)
    column_1 = args.column_1 - 1  ## Because Python index starts from 0
    column_2 = args.column_2 - 1
    outputfile = args.outputfile
    matchType = args.matchType

    if matchType == "partial":
        df = splitByQuantile_partial(file1, file2, column_1, column_2)
    elif matchType == "full":
        df = splitByQuantile_full(file1, file2, column_1, column_2)

    writeToFile(df, outputfile)
