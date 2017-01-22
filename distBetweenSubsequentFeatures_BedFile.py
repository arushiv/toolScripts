# Script function:
# 1. Given two bed files, output name_file_1, length_file_1, name_file_2, length_file_2, overlap 

#!/usr/bin/env python
import argparse
import csv
import pandas
import pybedtools

def replaceNegatives(x):
    if x <= 0:
        return 0
    else:
        return x


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="""Given a bed file, sort that file and compute the distance between each pair of subsequent features (matching chromosomes only). 
Usage:  python ~arushiv/toolScripts distBetweenSubsequentFeatures_BedFile.py a.bed""")
    parser.add_argument('bedfile', type=str, help="""bed file.""")
    parser.add_argument('outputfile', type=str, help="""Output file name.""")
    parser.add_argument('--full', action='store_const', const='full', help="""Output original bed file with an additional column at the end. Last feature of each chromosome is retained with 'NA' in the appended column. IF this option not specified, the output file will only contain 1 column with Distances, the last feature for each chromosome with NA in Distances will be removed.""")
    
    args = parser.parse_args()
     
    bedfile = pybedtools.BedTool(args.bedfile).sort().to_dataframe()

    if args.full != "full":       # Return only length column
        Distances = bedfile.groupby('chrom').apply(lambda x: x['start'].shift(-1) - x['end']).apply(lambda x: replaceNegatives(x)).dropna().astype(int)
        Distances.to_csv(args.outputfile, index=False)

    else:                         # Return full dataframe:
        bedfile.loc[:,'Distances'] = bedfile.groupby('chrom').apply(lambda x: x['start'].shift(-1) - x['end']).apply(lambda x: replaceNegatives(x)).reset_index().loc[:,0]
        # bedfile = bedfile.dropna()
        bedfile.loc[:,'Distances'] = bedfile.loc[:,'Distances'].astype(int, raise_on_error=False)
        bedfile.to_csv(args.outputfile, index=False, sep='\t', na_rep="na", header=False)
    
    
