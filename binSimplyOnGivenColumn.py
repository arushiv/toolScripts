#!/usr/bin/env python
import argparse
import csv
import glob
import os
import pandas



def filterPval(df, threshold):
    outdf = df[ (df['BetaAdjustedPval'] <= threshold) & (df['GeneType'] == "protein_coding") ]
    return  outdf


def absoluteBeta(df, column):
    df.loc[:,column] = abs(df[column])
    return df

def cutByQuantile(df, bins, column):
    # df['bin'] = pandas.qcut(df["absoluteBeta"], bins, labels=list(range(1, bins+1)))
    df.loc[:,'binNumber'] = pandas.qcut(df[column], bins, labels=list(range(1, bins+1))) #labels=list(range(1, (bins+1))))
    # df.loc[:,'bin'] = pandas.qcut(df[column], bins).str.replace(" ","-") #labels=list(range(1, (bins+1))))
    return df

def printdf(df, outputFile):
    df.to_csv(outputFile, sep='\t', index=False)

    
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="""Bin Dataframe based on one column, bin number appended in one columns. IMP: Script designed to filter for p value and take protein coding genes. Edit if this aspect not needed.)
Usage:  ~/myEnv/bin/python ~arushiv/toolScripts/binSimplyOnGivenColumn.py filename -b 4 outfile.txt """)
    parser.add_argument('dataframe', type=str, help="""Tab separated dataframe, with header""")
    parser.add_argument('-b', '--bins', type=int, help="""Number of bins""")
    parser.add_argument('-c', '--column', type=str, default='Slope', help="""Column name to bin on. Default = 'Slope'""")
    parser.add_argument('-abs', '--absolute', action="store_true", help="Bin after calculating absolute value of given column")
    parser.add_argument('outputFile', type=str, help="""Output file name""")
    parser.add_argument('-t', '--p_value_threshold', type=float, nargs='?', default=0.01, help="p value threshold to filter eqtl. Default = 0.01")
    parser.add_argument('-nf', '--noFilter', action="store_true", help="Override any pvalue or genetype filtering")
    parser.add_argument('-s', '--saveDifferent', action="store_true", help="save each bin into different file")

    args = parser.parse_args()
    threshold = args.p_value_threshold
    bins = args.bins
    outputFile = args.outputFile
    column = args.column
    
    df = pandas.read_csv(args.dataframe, sep='\t')

    if args.absolute:
        df = absoluteBeta(df, column)

    if args.noFilter:
        newdf = cutByQuantile(df, bins, column)
    else: 
        newdf = cutByQuantile(filterPval(df, threshold), bins, column)

    if args.saveDifferent:
        newdf.drop(column, axis=1, inplace=True)
        
        grouped = newdf.groupby(newdf.binNumber)
        for name, group in grouped:
            filename = outputFile + str(name)
            printdf(group, filename)
    else:
        printdf(newdf, outputFile)


