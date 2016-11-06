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

def subsetJasparJolma(bedfile):
    bedfile.columns = ['chrom','start','end','name','score','strand']
    jasparJolmaTFs = bedfile[bedfile['name'].str.contains("_") == False]
    encodeTFs = bedfile[bedfile['name'].str.contains("_")]
    return jasparJolmaTFs, encodeTFs

def splitByTF(encodeTFs):
    encodeTFs['motif'] = encodeTFs.name.str.split('_',1).str.get(1)              
    encodeTFs['name'] = encodeTFs.name.str.split('_',1).str.get(0)

    return encodeTFs

def bestScoringFeature(x):
    bedTf = pybedtools.BedTool.from_dataframe(x)                # Columns in input dataframe 'chr,  start,  end, name, score, strand,  motif'
    outdf = bedTf.intersect(bedTf.merge(), wao=True).to_dataframe()  # Columns in created dataframe are 'chrom  start  end name  score  strand thickStart  thickEnd  itemRgb blockCount blockSizes'
    # outdf['name'] = outdf.name.str.cat(outdf.thickStart.astype(str),sep='_')
    outdf['thickEnd'] = outdf['thickEnd'] + ":" + outdf['itemRgb'].map(str) + ":" + outdf['blockCount'].map(str) # thickEnd becomed 'key' to find max scoring motif
    outdf.drop(['itemRgb','blockCount','blockSizes'],1, inplace=True)
    return outdf[ outdf.groupby('thickEnd')['score'].transform(max) == outdf['score'] ].drop('thickEnd',1)

def printOutput(x, outputfile):
    x.to_csv(outputfile, header=False, sep='\t', index=False)

def operate(bedfile):
    jasparJolmaTFs, encodeTFs = subsetJasparJolma(bedfile)
    encodeTFs = splitByTF(encodeTFs)
    grouped = encodeTFs.groupby('name')
    # print grouped.apply(lambda x: bestScoringFeature(x))
    outdf = grouped.apply(lambda x: bestScoringFeature(x))

    outdf['name'] = outdf['name'] + "_" + outdf['thickStart'].map(str) 
    output = outdf.drop(['thickStart'],1).append(jasparJolmaTFs)
    output['score'] = output['score'].round().map(int)
    printOutput(pybedtools.BedTool.from_dataframe(output).sort().to_dataframe(), outputfile)

    
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="""If multiple motifs of one TF overlap each other (min 1 bp overlap), retain only the one with max score.\nFiltering performed only on Encode TFs Which contain '_'. String preceding first '_' is considered a TF.\n Requires pandas and pybedtools from python.\n Example ~/myEnv/bin/python removeOverlappingMotifs.py in.bed out.bed""")
    parser.add_argument('bedfile', type=str, help="""TF bed file, tab separated, should contain 6 columns: Chr,Start,end,name,score,strand. No header.""")
    parser.add_argument('outputfile', type=str, help="""Output file name. Will contain 4 columns: Chr;Start;end;'name,strand'.""")
    args = parser.parse_args()

    bedfile = pandas.read_csv(args.bedfile, header=None, sep='\t')
    outputfile = args.outputfile
    
    pandas.options.mode.chained_assignment = None

    operate(bedfile)
