# Script function: If multiple motifs of one TF overlap each other (min 1 bp overlap), retain only the one with max score.
# 1. Filtering performed only on Encode TFs Which contain '_'. String preceding first '_' is considered a TF.

#!/usr/bin/env python
import argparse
import csv
import gzip
import re
import operator
import glob
import math
import sys
import subprocess as sp
import os
import itertools
import shutil
import fnmatch
import pandas
import pybedtools

def subsetJaspar(bedfile):
    bedfile.columns = ['chrom','start','end','name','score']
    jasparTFs = bedfile.loc[bedfile.name.str.startswith('MA0')==True]
    encodeTFs = bedfile.loc[bedfile.name.str.startswith('MA0')==False]
    return jasparTFs, encodeTFs

def splitByTF(encodeTFs):
    newdf = encodeTFs['name'].apply(lambda x: pandas.Series(x.split('_',1)))
    newdf.columns = ['name','motif']
    formatdf = pandas.concat([encodeTFs[['chrom','start','end','score']], newdf], axis=1)
    return formatdf[['chrom','start','end','name','score','motif']]

def bestScoringFeature(x):
    bedTf = pybedtools.BedTool.from_dataframe(x)                # Columns in input dataframe 'chr,  start,  end, name, score, motif'
    mergedTf = bedTf.merge()                                     
    outdf = bedTf.intersect(mergedTf, wao=True).to_dataframe()  # Columns in created dataframe are 'chrom  start  end name  score  strand thickStart  thickEnd  itemRgb blockCount' 
    outdf['name'] = outdf.apply(lambda x:'%s_%s' % (x['name'],x['strand']),axis=1) 
    outdf['key'] = outdf.apply(lambda x:'%s:%s:%s' % (x['thickStart'], x['thickEnd'], x['itemRgb']), axis=1)
    return outdf[ outdf[['chrom', 'start', 'end', 'name', 'score', 'key']].groupby('key')['score'].transform(max) == outdf['score'] ][['chrom','start','end','name','score']]

def printOutput(x, outputfile):
    x.to_csv(outputfile, header=False, sep='\t', index=False)

def operate(bedfile):
    jasparTFs, encodeTFs = subsetJaspar(bedfile)
    formatdf = splitByTF(encodeTFs)
    grouped = formatdf.groupby('name')
    output = grouped.apply(lambda x: bestScoringFeature(x)).append(jasparTFs)
    printOutput(pybedtools.BedTool.from_dataframe(output).sort().to_dataframe(), outputfile)

    
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="""If multiple motifs of one TF overlap each other (min 1 bp overlap), retain only the one with max score.\nFiltering performed only on Encode TFs Which contain '_'. String preceding first '_' is considered a TF.\n Requires pandas and pybedtools from python.\n Example ~/myEnv/bin/python removeOverlappingMotifs.py in.bed out.bed""")
    parser.add_argument('bedfile', type=str, help="""TF bed file, tab separated, should contain 5 columns: Chr,Start,end,name,score. No header.""")
    parser.add_argument('outputfile', type=str, help="""Output file name.""")
    args = parser.parse_args()

    bedfile = pandas.read_csv(args.bedfile, header=None, sep='\t')
    outputfile = args.outputfile
    
    operate(bedfile)
