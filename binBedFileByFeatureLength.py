#!/usr/bin/env python
import argparse
import csv
import glob
import os
import pandas
import ggplot as gp  
import pybedtools as pb

def printFunc(df, outputFile):
        name = "%s.%s.bed" %(outputFile, df.loc[0]['bin'] + 1)
        df.loc[:,['chrom','start','end','bin']].to_csv(name, sep="\t", index=False, header=False)
    

def cutByLength(bedfile, bins, outputFile):
    bedfile.loc[:,'lengthFeature'] = bedfile.loc[:,'end'] - bedfile.loc[:,'start']
    bedfile.loc[:,'bin'] = pandas.qcut(bedfile.loc[:,'lengthFeature'], bins, labels=False)

    bedfile.groupby('bin').apply(lambda x: printFunc(x.reset_index(), outputFile))


    
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="""Bin a bed file based on length of features.)
Usage:  python ~arushiv/toolScripts/binBedFileByFeatureLength.py <bedfile> <outname> -b <bins>""")
    parser.add_argument('bedfile', type=str, help="""Bed file to be split into bins by length of features""")
    parser.add_argument('-b', '--bins', type=int, help="""Number of bins""")
    parser.add_argument('outputFile', type=str, help="""Output file identifier, each binned file will be named as outputFile.binNumber.bed""")


    args = parser.parse_args()
    bins = args.bins
    outputFile = args.outputFile
    
    bedfile = pb.BedTool(args.bedfile).sort().merge().to_dataframe()

    cutByLength(bedfile, bins, outputFile)
