# Script function:
# 1. From the GREGOR output folder, obtain GWAS loci that overlap. Take LD SNPs for each  GWAS locus .
# 2. Overlap all tag + LD with the original bed file, report all overlapping SNPs


#!/usr/bin/env python
import argparse
import glob
import os
import pandas
import numpy as np
from StringIO import StringIO
import pybedtools
                        
def getLdSnps(files, gregorOutputFolder):  # from index.SNP.in.LD.Chrxx files, make dataframe of overlapping index snps and their LD buddies
        name = os.path.basename(files.rstrip())
        mdf = pandas.DataFrame()
        for overlapFile in glob.glob(os.path.join(gregorOutputFolder, name, "index.snp.LD.on.chr*")):
                # if sum(1 for line in open(overlapFile)) > 1:
                df = pandas.read_table(overlapFile, sep='\t', dtype={'LD_buddy_pos': str}) 
                mdf = mdf.append(df, ignore_index=True)
                mdf['filename'] = name
        return mdf
        

def makeLdBedFile(df):  # From index SNP and LD buddy dataframe, make bed file with chrom start end for each snp and LD snp

        s = df['LD_buddy_pos'].str.split("|", expand=True).apply(pandas.Series, 1).stack(dropna=True)
        s.index = s.index.droplevel(-1)
        s.name = 'pos'
        dreturn = df.join(s)
        dreturn["chrom"], dreturn["indexSnp"] = zip(*dreturn["index_SNP"].str.split(':').tolist())
        del dreturn["index_SNP"]
        dreturn = dreturn[["chrom","pos","pos","indexSnp","filename"]]
        dreturn.iloc[:,1] = dreturn.iloc[:,1].apply(int) - 1

        return dreturn


    
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='From enrichment in features GREGOR output, retrieve LD buddies of index_SNPs that overlap. Working directory to be the one which contains bed_path_file if it is provided in an argument', usage='python gregorFindAllLdSnps.py output_T2D output.txt -p bedfile.txt')
    parser.add_argument('gregorOutputFolder', type=str, help="""GREGOR output directory.""")
    parser.add_argument('-p','--bed_path_file', type=str, help="""The bed path file used for running GREGOR.""")
    parser.add_argument('-f','--feature', type=str, help="""Only fetch overlapping SNPs for a specific bed file paths. Paths supplied as comma separated strings""")
    parser.add_argument('outputfile', help ="""Output file. Will contain these new columns: chrom_overlapSNP start_overlapSNP end_overlapSNP pos_indexSNP chrom_bedStart start_bedStart end_bedStart filename""")
    args = parser.parse_args()

    gregorOutputFolder = args.gregorOutputFolder
    outputfile = args.outputfile
    outdf = pandas.DataFrame()

    if args.bed_path_file is not None:
        with open(args.bed_path_file) as f:
                bedList = f.readlines()
    elif args.feature is not None:
            bedList = args.feature.split(",")
    else:
            print "Please provide at least one argument -p or -f"


for files in bedList:
        df = getLdSnps(files, gregorOutputFolder)
        if not df.empty:
                print files
                ldBedFile = makeLdBedFile(df)
                outdf = outdf.append(ldBedFile, ignore_index=True)

outdf.to_csv(outputfile, sep='\t', index=False, header=False)
