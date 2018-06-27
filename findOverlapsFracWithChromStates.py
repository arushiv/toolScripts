

#!/usr/bin/env python
import argparse
import csv
import glob
import os
import pandas
import pybedtools


def lengthOfBedFile(x, lengthFile1):
    return "{0:.2f}".format(float((x['end'] - x['start']).sum(axis=0))/lengthFile1)

def getname(s):
    return os.path.splitext(os.path.basename(s))[0].replace('.bed','')

def overlapFrac(x, bedfile2):    # Intersect each row of bedfile 1 with bedfile 2, output fraction of overlap with each unique 'name' column value of bedfile 2
    intersect = bedfile2.intersect(pybedtools.BedTool.from_dataframe(x)).to_dataframe()
    lengthFile1 = x.loc[0]['end'] - x.loc[0]['start']
    df = intersect.groupby('name').apply(lambda x: lengthOfBedFile(x,lengthFile1)).to_frame().reset_index()
    return ';'.join((df['name'] + ',' + df[0].map(str)).tolist())


def computeOverlap(bedfile1, bedfile2name):
    bedfile2 = pybedtools.BedTool(bedfile2name)
    nameFile = getname(bedfile2name)
    for index, row in bedfile1.iterrows():
        bedfile1.loc[index, nameFile] = overlapFrac(row.to_frame().transpose().reset_index().drop('index',axis=1), bedfile2)
    return bedfile1


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="""For each row in bedfile 1, find fraction of overlaps with bedfile 2 grouped by unique name field in bedfile 2. Output file will contain one extra column over bedfile 1 noting the fraction of overlap with each chromatin state separated by semicolons.
Usage: python ~arushiv/toolScripts/findOverlapsFracWithChromStates.py file1.bed file2.bed out.bed""")
    parser.add_argument('bedfile1', type=str, help="""bed file 1. Each line of this bed file will be tested for overlaps. IMP: Supply header, first 3 columns should be chrom, start, end """)
    parser.add_argument('--header', nargs='+', help="""If bed file 1 does not contain header, provide a space separated list. IMP first 3 columns should be chrom, start, end """)
    # parser.add_argument('-d', '--directory', type=str, default='/lab/work/arushiv/13ChromatinStatesByCellType', help="""Take all bed.gz files from this directory.""")
    parser.add_argument('bedfile2', type=str, help="""bed file 2.""")
    parser.add_argument('outputfile', help="""Output file name.""")

    args = parser.parse_args()

    if args.header is not None:
        bedfile1 = pandas.read_csv(args.bedfile1, sep='\t', header=None, names=args.header)
    else:
        bedfile1 = pandas.read_csv(args.bedfile1, sep='\t')
        
    bedfile2name = args.bedfile2

    computeOverlap(bedfile1, bedfile2name).to_csv(args.outputfile, sep='\t', index=False)
