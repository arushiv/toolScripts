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
    # ndf = df.closest(pybedtools.BedTool(bedfile), d=True).to_dataframe() # closest for gene
    ndf = pybedtools.BedTool(bedfile).closest(df, d=True).to_dataframe()  # closest for annotation
    ndf = ndf[ndf.columns[-1:]]
    ndf = ndf[ndf.iloc[:,0] != -1]
    ndf.loc[:,'fileName'] = getname(bedfile)
    return ndf


def getNearestDistance(df, bedfilelist):
    odf = pandas.concat([fixFunc(df, bedfile) for bedfile in bedfilelist])
    return odf
    # [bedfile.closest(df, d=True).to_dataframe() for bedfile in bedfilelist]

def getname(s):
    return os.path.splitext(os.path.basename(s))[0].split('.')[1]


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="""1. Given a  bed files, randmly sample on each unique item in name or other column. 2. Go back to original bedfile and fetch all coordinates corresponding to sampled name field. 3. Compute bedtools closest given a list of annotation files for each sampling iteration and return dataframe. Usage:  ~/myEnv/bin/python ~arushiv/toolScripts/randomlySampleGenesAndComputeNearestDistance.py  TssFile.bed <list.of.annotation.BedFiles> -n <number of sampling iterations>""")
    parser.add_argument('tssFile', type=str, help="""bed file on which sampling has to be done. No header""")
    parser.add_argument('bedfilelist', nargs='+', type=str, help="""List of annotation bedfiles [3 column tab separated no header] to compute closest distance after each sampling""")
    parser.add_argument('-o', '--outfile', type=str, default="outfile.dat", help="""Output file name.""")
    args = parser.parse_args()

    
    tssFile = pybedtools.BedTool(args.tssFile).sort()
    bedfilelist = args.bedfilelist
   
    outdf = getNearestDistance(tssFile, bedfilelist)
    outdf.to_csv(args.outfile, index=False, sep='\t')
