# Script function:
# 1. Given a  bed files, randmly sample on each unique item in name or other column. 2. Go back to original bedfile and fetch all coordinates corresponding to sampled name field. 3. Compute bedtools closest given a list of annotation files for each sampling iteration and return dataframe.


#!/usr/bin/env python
import argparse
import glob
import os
import pandas
import pybedtools
import numpy
from statsmodels.distributions.empirical_distribution import ECDF


def fixFunc(df, bedfile):
    # ndf = df.closest(pybedtools.BedTool(bedfile), d=True).to_dataframe()  # Closest for genes
    ndf = pybedtools.BedTool(bedfile).closest(df, d=True).to_dataframe()  # Closest for annotations
    ndf = ndf[ndf.columns[-1:]]  # Retain only last column of the closest distance
    ndf.columns = ['distance']
    ndf = ndf[ndf.iloc[:,0] != -1]  # pybedtools.closest returns -1 if no feature is found on the same chromosome. Remove this from the output dataframe
    namelist = getname(bedfile)
    ndf.loc[:,'cell'], ndf.loc[:,'annotation'] = namelist[0], namelist[1]
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
    return os.path.splitext(os.path.basename(s))[0].split('.')

def ecdfOfDataframeColumn(x):
    ecdf = ECDF(x)
    rangeOfx = numpy.arange(0, 8.5, step=0.05)
    return ecdf(rangeOfx)

def getXColumn(x):
    x.loc[:,'xvals'] = numpy.arange(0, 8.5, step=0.05).tolist()
    return x
    
def getEcdf(df, groupings):
    df.loc[:,'distance'] = numpy.log10(df['distance'] + 1)
    fulldf = df.groupby(groupings)['distance'].apply(pandas.core.groupby.Series.tolist).apply(ecdfOfDataframeColumn).reset_index()
    fulldf = fulldf.rename(columns={'distance': 'ecdf_y'})
    fulldf = fulldf['ecdf_y'].apply(pandas.Series).stack().reset_index(level=1, drop=True).to_frame('ecdf_y').join(fulldf[groupings], how='left')
    fulldf = fulldf.groupby(groupings).apply(getXColumn)
    return fulldf


def getMean(df, groupings):
    outdf = df.groupby(groupings).agg([numpy.mean, numpy.std]).reset_index()
    outdf.columns = outdf.columns.droplevel()
    groupings.extend(['ecdf_y_shuffleMean','ecdf_y_shuffleStd'])
    outdf.columns = groupings
    return outdf

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="""Compute enrichment for nearest distance to TSS with genes by quintile with mean of random shuffles done n times. Shuffling: Given a  bed file, randomly sample on each unique item in name or other column. 2. Go back to original bedfile and fetch all coordinates corresponding to sampled name field. 3. Compute bedtools closest given a list of annotation files for each sampling iteration and return dataframe. Usage:  ~/myEnv/bin/python ~arushiv/toolScripts/randomlySampleGenesAndComputeNearestDistance.py  TssFile.bed <list.of.annotation.BedFiles> -n <number of sampling iterations>""")
    parser.add_argument('tssFile', type=str, help="""bed file on which sampling has to be done. No header""")
    parser.add_argument('bedfilelist', nargs='+', type=str, help="""List of annotation bedfiles [3 column tab separated no header] to compute closest distance after each sampling""")
    parser.add_argument('-df','--distanceAnnotFile', default='/lab/work/arushiv/chromatin/regulatory_regions_review/distanceDistFromTSS/nearestDistFromGencodeV19TSSUniq_geneTssByLclEsi.dat.gz', help="""Data file with nearest gene TSS distance values for each regulatory annotation. Defualt = '/lab/work/arushiv/chromatin/regulatory_regions_review/distanceDistFromTSS/nearestDistFromGencodeV19TSSUniq_geneTssByLclEsi.dat.gz'""")
    parser.add_argument('-s', '--numberOfGenes', type=int, default=3700, help="""Number of genes in each sampling. i.e. number of genes in each quintile. Default 3700.""")
    parser.add_argument('-c', '--colNumber', type=int, default=4, help="""Column number of gene names from which unique genes have to be sampled""")
    parser.add_argument('-r', '--runs', type=int, default=5, help="""Number of sampling runs to perform. Data will be recorded with additional column of the run number. Default = 5 runs""")
    parser.add_argument('-o', '--outfile', type=str, default="outfileshuffle.dat", help="""Output file name.""")
    args = parser.parse_args()

    
    tssFile = pandas.read_csv(args.tssFile, header=None, sep='\t')
    bedfilelist = args.bedfilelist
    numberOfGenes = args.numberOfGenes
    colNumber = args.colNumber - 1   # python index from 0

    distanceAnnotFile = pandas.read_csv(args.distanceAnnotFile, sep='\t', header=0, names=['cell','annotation','quintile','distance'])
    groupings_distAnnot = ['cell','annotation','quintile']
    distanceAnnotFile_ecdf = getEcdf(distanceAnnotFile, groupings_distAnnot)   #.to_csv("test.df", sep='\t', index=False)
    
    # Shufflings
    output = []
    for runNumber in range(args.runs):
        print runNumber + 1
        tdf = sampleUnique(tssFile, colNumber, numberOfGenes)
        distDf = getNearestDistance(tdf, bedfilelist, runNumber)
        output.append(distDf)

    shuffledf = pandas.concat(output)

    groupings_shuffle = ['cell','annotation','runNumber']
    groupings_shuffleMeans = ['cell','annotation','xvals']
    shuffle_MeanEcdf = getMean(getEcdf(shuffledf, groupings_shuffle), groupings_shuffleMeans)

    pandas.merge(distanceAnnotFile_ecdf, shuffle_MeanEcdf, on=['cell','annotation','xvals']).to_csv(args.outfile, index=False, sep='\t', header=True)

    
    # outputdf.to_csv(args.outfile, index=False, sep='\t', header=True)

    # compile data into bins
#### Function to bin data and record mean and sd #####
# def binData(df):
#     df.loc[:,'distance'] = numpy.log10(df['distance'] + 1)
#     rangeOfx = numpy.arange(0, 8, step=0.05).tolist()
#     labels = list((x + 0.025 for x in rangeOfx[:-1]))
#     grouped = df.groupby(['cell','annotation','runNumber'])
#     outdf = []
#     for name, group in grouped:
#         group.loc[:,'binNumber'] = pandas.cut(group['distance'], bins=rangeOfx, labels=labels)
#         ndf = group.groupby(['cell','annotation','runNumber','binNumber']).agg([numpy.mean, numpy.std]).reset_index()
#         outdf = outdf.append(ndf)
#     outdf = pandas.concat(outdf)
#     return outdf
   



