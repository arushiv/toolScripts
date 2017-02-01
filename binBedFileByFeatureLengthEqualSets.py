#!/usr/bin/env python
import argparse
import csv
import glob
import os
import pandas
import pybedtools as pb

def printFunc(df, outputFile):
        # name = "%s.%s.bed" %(outputFile, df.loc[0]['bin'] + 1)
        # df.loc[:,['chrom','start','end']].to_csv(name, sep="\t", index=False, header=False)
        df.to_csv(outputFile, sep='\t', index=False, header=False)
def func(x, bins):
        if x > bins:
                return x-1
        else:
                return x

def cutByLength(bedfile, bins):
    bedfile[bedfile.groupby(['chrom','start','end'])['score'].transform(max) == bedfile['score']]    # For features with same coordinates, filter for one with max score
    bedfile.loc[:,'lengthFeature'] = bedfile.loc[:,'end'] - bedfile.loc[:,'start']                   # Calculate length of feature
    bedfile.sort_values(by='lengthFeature', inplace=True)                                            # Sort by length
    estimate = bedfile['lengthFeature'].sum()/bins                        # Approx. size of each bin
    bedfile.loc[:,'cumsum'] = bedfile.loc[:,'lengthFeature'].cumsum()     # Cumulative sum to calculate length of each bin 
    bedfile.loc[:,'bin'] = bedfile.loc[:,'cumsum'].apply(lambda x: int(x/(bedfile.iloc[-1]['cumsum']/bins))+1).apply(lambda x: func(x, bins))    ## Us cumulative sum to assign bins
    del bedfile['cumsum']
    return bedfile
    # bedfile.loc[:,'bin'] = pandas.qcut(bedfile.loc[:,'lengthFeature'], bins, labels=False)

    # d1 = bedfile.loc[:,['cumsum','bin']].groupby('bin').max()
    # print d1
    # print d1['cumsum'].shift(-1) - d1['cumsum'] 

    
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="""1. For a bed file, in case of same segments with differring scores, retain only segments that have the highest scores. 2. Sort bedfile by length of features. 3. Add a column that records bin number, keeping total length of features in each bin approx. equal. Usage:  python ~arushiv/toolScripts/binBedFileByFeatureLength.py <bedfile> <outname> -b <bins>""")
    parser.add_argument('bedfile', type=str, help="""Bed file to be split into bins by length of features""")
    parser.add_argument('-b', '--bins', type=int, help="""Number of bins""")
    parser.add_argument('outputFile', type=str, help="""Output file identifier, each binned file will be named as outputFile.binNumber.bed""")


    args = parser.parse_args()
    bins = args.bins
    outputFile = args.outputFile
    
    bedfile = pb.BedTool(args.bedfile).sort().to_dataframe()

    printFunc(cutByLength(bedfile, bins), outputFile)
