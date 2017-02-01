# Script function:
# 1. Bin the eqtls using a sliding window of beta values 
# Colnames: ['GeneName', 'Strand', 'GencodeLevel', 'GeneType', 'GeneID', 'ChrPheno', 'StartPheno', 'EndPheno', 'BestExonID', 'NumExons', 'NumVariantCis', 'DistanceWithBest', 'SNPid', 'A1', 'A2', 'MAF', 'SNPchr', 'StartSNP', 'EndSNP', 'Nominal_Pval', 'Slope', 'EmpiricalAdjustedPval', 'BetaAdjustedPval', 'eQTLnum', 'Backward_eQTLnum']
#!/usr/bin/env python
import argparse
import csv
import glob
import os
import pandas
import ggplot as gp  


def filterPval(df, threshold):
    outdf = df[ (df['BetaAdjustedPval'] <= threshold) & (df['GeneType'] == "protein_coding") ]
    return  outdf


def absoluteBeta(df):
    df.loc[:,'absoluteBeta'] = abs(df['Slope'])
    return df

def cutByQuantile(df, bins):
    # df['bin'] = pandas.qcut(df["absoluteBeta"], bins, labels=list(range(1, bins+1)))
    df.loc[:,'bin'] = pandas.qcut(df["MAF"], bins,).str.replace(" ","-") #labels=list(range(1, (bins+1))))
    return df

def printdf(df, outputFile):
    df.to_csv(outputFile, sep='\t', index=False)

def plotDensity(df):
    p = gp.ggplot(df, gp.aes(x='absoluteBeta', fill='bin')) +\
    gp.geom_density() +\
    scale_x_log10() +\
    gp.theme_bw()

    p.save("BetaByMaf.pdf")



    
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="""Bin Dataframe based on one column, bin number appended in one columns. IMP: Script designed to filter for p value and take protein coding genes. Edit if this aspect not needed.)
Usage:  ~/myEnv/bin/python ~arushiv/toolScripts/binSimplyOnGivenColumn.py filename -b 4 outfile.txt """)
    parser.add_argument('dataframe', type=str, help="""Tab separated dataframe, with header""")
    parser.add_argument('-b', '--bins', type=int, help="""Number of bins""")
    parser.add_argument('outputFile', type=str, help="""Output file name""")
    parser.add_argument('-t', '--p_value_threshold', type=float, nargs='?', default=0.01, help="p value threshold to filter eqtl. Default = 0.01")

    args = parser.parse_args()
    threshold = args.p_value_threshold
    bins = args.bins
    outputFile = args.outputFile
    
    df = pandas.read_csv(args.dataframe, sep='\t')

    newdf = cutByQuantile(absoluteBeta(filterPval(df, threshold)), bins)
    printdf(newdf, outputFile)
    # plotDensity(newdf)

