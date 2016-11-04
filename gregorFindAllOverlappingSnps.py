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
                df = pandas.read_table(overlapFile, sep='\t') 
                mdf = mdf.append(df, ignore_index=True)

        return mdf
        
def makeLdBedFile(df):  # From index SNP and LD buddy dataframe, make bed file with chrom start end for each snp and LD snp
        snp_df = pandas.DataFrame()
        for index, row in df.iterrows():
                chrom = row['index_SNP'].split(":")[0]
                pos = row['LD_buddy_pos'].split("|")
                for stuff in pos:
                        d = pandas.read_table(StringIO("%s\t%s\t%s"%(chrom,stuff,stuff)), sep="\t", header=None)
                        d['index_SNP'] = row['index_SNP']
                        # d = pandas.read_table(StringIO("%s\t%s\t%s"%(chrom,int(stuff)-1,stuff)), sep="\t", header=None)
                        snp_df = snp_df.append(d, ignore_index=True)
        return snp_df

def determineOverlaps(df, files): # Intersect original bed file and the snp bed file with snp and LD snps
        snp_bed = pybedtools.BedTool.from_dataframe(df)
        bedfile = pybedtools.BedTool(files.rstrip())
        output = snp_bed.intersect(bedfile).to_dataframe()

        return output

    
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='From enrichment in features GREGOR output, fetch all tag + LD SNP overlaps for a given bed file. Working directory to be the one which contains bed_path_file if it is provided in an argument', usage='~/myEnv/bin/python gregorFindAllOverlappingSnps.py <gregor output folder> <output file> -f <comma separated paths to bed files to check overlaps in>')
    parser.add_argument('gregorOutputFolder', type=str, help="""GREGOR output directory.""")
    parser.add_argument('-p','--bed_path_file', type=str, help="""The bed path file used for running GREGOR.""")
    parser.add_argument('-f','--feature', type=str, help="""Only fetch overlapping SNPs for a specific bed file paths. Paths supplied as comma separated strings""")
    parser.add_argument('outputfile', help ="""Output file. Will contain these new columns: chrom   pos index_SNP filename""")
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
                overlappingSnps = determineOverlaps(ldBedFile, files)
                overlappingSnps['filename'] = os.path.basename(files.rstrip())
                outdf = outdf.append(overlappingSnps[['chrom','end','name','filename']], ignore_index=True)
        else:
                continue

        
outdf.to_csv(outputfile, sep='\t', index=False, header=["snp_chrom","snp_pos","index_snp","filename"])


