# Script function:
# 1. From the GREGOR output folder, find enriched features based on input pvalue thresholds and/or input enrichment cutoff.
# 2. For qualifying features, find eQTL SNPs that overlap.
# 3. For overlapping SNPs, go to the dataframe and fetch gene list

#!/usr/bin/env python
import argparse
import glob
import os
import pandas
import numpy as np


def getStatsDataFrame(gregorOutputFolder):
    stats = pandas.read_csv(os.path.join(gregorOutputFolder, 'StatisticSummaryFile.txt'), sep='\t')
    return stats

    
def findSignificantFeatures(stats, pval, minEnrich):
    featureFrame = stats.loc[(stats['PValue'] <= pval) & (np.log2(stats['InBed_Index_SNP'] / stats['ExpectNum_of_InBed_SNP']) >= minEnrich)]
    return featureFrame

def fetchFrameForSpecific(stats, tfstring):
    featureFrame = pandas.DataFrame()
    for tf in tfstring.split(","):
        featureFrame = pandas.concat([featureFrame, stats[stats.Bed_File.str.contains(tf)]], ignore_index=True)
    return featureFrame

def getOverlappingSnps(featureFrame, gregorOutputFolder):
    df = pandas.DataFrame()
    for feature in featureFrame['Bed_File']:
        for filename in glob.glob(os.path.join(gregorOutputFolder, feature, "index.snp.LD.on.chr*")):
            readdf = pandas.read_csv(filename, sep='\t', usecols=['index_SNP'])
            readdf['Bed_File'] = feature
            df = pandas.concat([df, readdf], ignore_index=True)
    return df[['Bed_File','index_SNP']]
        
def fetchGenesForOverlappingSnps(df, eqtl_gene_file):
    return pandas.merge(df, eqtl_gene_file, how='inner', on=['index_SNP']).sort_values(by='Bed_File')

def printOutput(x, outputfile):
    x.to_csv(outputfile, sep='\t', index=False)
    
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='From eQTL enrichment in features GREGOR output, fetch genes associated with eQTL')
    parser.add_argument('gregorOutputFolder', type=str, help="""GREGOR output directory.""")
    parser.add_argument('eqtl_gene_file', type=str, help="""The tab separated file with eQTL variants and associated genes. First 2 columns should be of the format Chr:Pos, eg. "chr16:81820708\t7SK", GeneName, with header names ['index_SNP','gene']. Extra columns accepted, will be appended in the output file.""")
    parser.add_argument('-p','--pvalThreshold', type=float, nargs='?', default=0.05,
                        help="""pvalue threshold; default 0.05""")
    parser.add_argument('-e','--enrichmentThreshold', type=float, nargs='?', default=0.0, help="""Minimum log2(fold enrichment) enrichment threshold. Default 0.0""")
    parser.add_argument('-f','--feature', type=str, default='none', help="""Submit specific features to fetch genes for, as comma separated strings. A string find will be performed on folders in GREGOR output to fetch overlaps""")
    parser.add_argument('outputfile', help ="""Output file. Will contain atleast: Bed_File,index_SNP,gene columns""")
    args = parser.parse_args()

    gregorOutputFolder = args.gregorOutputFolder
    
    eqtl_gene_file = pandas.read_csv(args.eqtl_gene_file, sep='\t')#, header=None, names=['index_SNP','gene'])
    outputfile = args.outputfile
    
    pval = args.pvalThreshold
    minEnrich = args.enrichmentThreshold
    tfstring = args.feature
        
    stats = getStatsDataFrame(gregorOutputFolder)
    if tfstring == "none":
        featureFrame = findSignificantFeatures(stats, pval, minEnrich)
    else:
        featureFrame = fetchFrameForSpecific(stats, tfstring)
        
    overlappingSNPs = getOverlappingSnps(featureFrame, gregorOutputFolder)
    printOutput(fetchGenesForOverlappingSnps(overlappingSNPs, eqtl_gene_file), outputfile)
