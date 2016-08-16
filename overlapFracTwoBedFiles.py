# Script function:
# 1. Given two bed files, output name_file_1, length_file_1, name_file_2, length_file_2, overlap 

#!/usr/bin/env python
import argparse
import csv
import glob
import os
import pandas
import pybedtools


def lengthOfBedFile(x):
    if x.count() == 0:
        return 0
    else:    
        df = x.to_dataframe()
        return (df['end'] - df['start']).sum(axis=0)

def getname(s):
    return os.path.splitext(os.path.basename(s))[0]


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="""Given two bed files, output name_file_1, length_file_1, name_file_2, length_file_2, overlap.
Usage:  ~/myEnv/bin/python ~arushiv/toolScripts/overlapFracTwoBedFiles.py a.bed b.bed""")
    parser.add_argument('bedfile1', type=str, help="""bed file 1.""")
    parser.add_argument('bedfile2', type=str, help="""bed file 2.""")
    # parser.add_argument('outputfile', type=str, help="""Output file name. Will contain 4 columns: Chr;Start;end;'name,strand'.""")
    args = parser.parse_args()

    
    
    bedfile1 = pybedtools.BedTool(args.bedfile1)
    bedfile2 = pybedtools.BedTool(args.bedfile2)
    intersect = bedfile1.intersect(bedfile2)

    print getname(args.bedfile1), lengthOfBedFile(bedfile1), getname(args.bedfile2), lengthOfBedFile(bedfile2), lengthOfBedFile(intersect)
    
    
