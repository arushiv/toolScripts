# Script function:
# 1. Given a  bed files, randmly sample on each unique item in name or other column. 2. Go back to original bedfile and fetch all coordinates corresponding to sampled name field. 3. Compute bedtools closest given a list of annotation files for each sampling iteration and return dataframe.


#!/usr/bin/env python
import argparse
import csv
import glob
import os
import pandas
import pybedtools


def fixFunc(df, bedfile):
    # ndf = df.closest(pybedtools.BedTool(bedfile), d=True).to_dataframe()  # Closest for genes
    ndf = pybedtools.BedTool(bedfile).closest(df, d=True).to_dataframe()  # Closest for annotations
    ndf = ndf[ndf.columns[-1:]]  # Retain only last column of the closest distance
    ndf = ndf[ndf.iloc[:,0] != -1]  # pybedtools.closest returns -1 if no feature is found on the same chromosome. Remove this from the output dataframe
    ndf.loc[:,'fileName'] = getname(bedfile)
    return ndf

def sampleUnique(tssFile, colNumber, numberOfGenes):
    genelist = tssFile.iloc[:,colNumber].drop_duplicates().sample(numberOfGenes)
    return pybedtools.BedTool.from_dataframe(tssFile[tssFile.iloc[:,colNumber].isin(genelist)]).sort()

def getNearestDistance(df, bedfilelist, runNumber):
    returndf = pandas.concat([fixFunc(df, bedfile) for bedfile in bedfilelist])
    returndf.loc[:,'runNumber'] = "sh" + str(runNumber + 1)
    return returndf
    # [bedfile.closest(df, d=True).to_dataframe() for bedfile in bedfilelist]

def getname(s):
    return os.path.splitext(os.path.basename(s))[0].split('.')[1]


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="""1. Given a  bed files, randmly sample on each unique item in name or other column. 2. Go back to original bedfile and fetch all coordinates corresponding to sampled name field. 3. Compute bedtools closest given a list of annotation files for each sampling iteration and return dataframe. Usage:  ~/myEnv/bin/python ~arushiv/toolScripts/randomlySampleGenesAndComputeNearestDistance.py  TssFile.bed <list.of.annotation.BedFiles> -n <number of sampling iterations>""")
    parser.add_argument('tssFile', type=str, help="""bed file on which sampling has to be done. No header""")
    parser.add_argument('bedfilelist', nargs='+', type=str, help="""List of annotation bedfiles [3 column tab separated no header] to compute closest distance after each sampling""")
    parser.add_argument('-s', '--numberOfGenes', type=int, default=3700, help="""Number of genes in each sampling. Default 3700.""")
    parser.add_argument('-c', '--colNumber', type=int, default=4, help="""Column number of gene names from which unique genes have to be sampled""")
    parser.add_argument('-r', '--runs', type=int, default=5, help="""Number of sampling runs to perform. Data will be recorded with additional column of the run number. Default = 5 runs""")
    parser.add_argument('-o', '--outfile', type=str, default="outfileshuffle.dat", help="""Output file name.""")
    args = parser.parse_args()

    
    tssFile = pandas.read_csv(args.tssFile, header=None, sep='\t')
    bedfilelist = args.bedfilelist
    numberOfGenes = args.numberOfGenes
    colNumber = args.colNumber - 1   # python index from 0

    outputdf = pandas.DataFrame()
    
    for runNumber in range(args.runs):
        print runNumber + 1
        tdf = sampleUnique(tssFile, colNumber, numberOfGenes)
        outdf = getNearestDistance(tdf, bedfilelist, runNumber)
        outputdf = outputdf.append(outdf, ignore_index=True)

    outputdf.to_csv(args.outfile, index=False, sep='\t', header=True)
