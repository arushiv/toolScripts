# Script function:
# 1. Bin the eqtls using a sliding window of beta values 
# Colnames: ['GeneName', 'Strand', 'GencodeLevel', 'GeneType', 'GeneID', 'ChrPheno', 'StartPheno', 'EndPheno', 'BestExonID', 'NumExons', 'NumVariantCis', 'DistanceWithBest', 'SNPid', 'A1', 'A2', 'MAF', 'SNPchr', 'StartSNP', 'EndSNP', 'Nominal_Pval', 'Slope', 'EmpiricalAdjustedPval', 'BetaAdjustedPval', 'eQTLnum', 'Backward_eQTLnum']
#!/usr/bin/env python
import argparse
import csv
import glob
import os
import pandas
# import pybedtools


def filterPval(df, threshold):
    outdf = df[ (df['BetaAdjustedPval'] <= threshold) & (df['GeneType'] == "protein_coding") ]
    return  outdf


def absoluteBeta(df):
    df.loc[:,'absoluteBeta'] = abs(df['Slope'])
    return df

def retainUniqueSnpsByPval(df):
    outdf = df.sort(['BetaAdjustedPval'], ascending = 1)
    outdf = outdf.drop_duplicates(subset=['SNPchr','EndSNP'], keep = 'first')
    print outdf
    return outdf

def cutByQuantile(df, bins):
    # df['bin'] = pandas.qcut(df["absoluteBeta"], bins, labels=list(range(1, bins+1)))
    df.loc[:,'slidingBin'] = pandas.qcut(df["absoluteBeta"], (2*bins), labels=list(range(1, (2*bins)+1)))
    return df

def printBySlidingBins(df, bins, identifier):
    slidingBins = 2*bins
    outDict = {}
    for slidingBin in list(range(1,slidingBins)):
        filename = identifier + "." + str(slidingBin) + ".dat"
        outdf = df[(df['slidingBin'] == slidingBin) | (df['slidingBin'] == slidingBin + 1) ]
        outdf.to_csv(filename, sep='\t', index=False)
        # outDict[slidingBin] = outdf

    # print outDict
        
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="""Bin eQTLs by value of slope or beta using a sliding window. Eg. If supplied bin parameter is 4, 7 (2*4 - 1) bins with sliding window equivalent to 8 bins will be created.)
Usage:  ~/myEnv/bin/python ~arushiv/toolScripts/bineQTLByBeta.py filename -b 4 -o folder/eqtl""")
    parser.add_argument('dataframe', type=str, help="""Tab separated dataframe, with header""")
    parser.add_argument('-b', '--bins', type=int, help="""Number of bins""")
    parser.add_argument('-o', '--outputFileIdentifier', type=str, nargs='?', default="eqtl", help="""Output file name identifier. Eg. identifier folder/eqtl will result in files: folder/eqtl.1.dat; folder/eqtl.2.dat and so on.""")
    parser.add_argument('-t', '--p_value_threshold', type=float, nargs='?', default=0.01, help="p value threshold to filter eqtl. Default = 0.01")

    args = parser.parse_args()
    threshold = args.p_value_threshold
    bins = args.bins
    identifier = args.outputFileIdentifier
    
    df = pandas.read_csv(args.dataframe, sep='\t')

    newdf = absoluteBeta(filterPval(df, threshold))
    retainUniqueSnpsByPval(newdf)
    # printBySlidingBins(cutByQuantile(newdf, bins), bins, identifier)
    
