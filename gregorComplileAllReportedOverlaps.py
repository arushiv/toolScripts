# Script function:
# 1. From the GREGOR output folder, obtain index_SNPs and/LD buddy that is reported to overlap. 


#!/usr/bin/env python
import argparse
import glob
import os
import pandas
import numpy as np
import pybedtools
                        
def getLdSnps(files, gregorOutputFolder):  # from index.SNP.in.LD.Chrxx files, make dataframe of overlapping index snps and their LD buddies
        name = os.path.basename(files.rstrip())
        mdf = pandas.DataFrame()
        for overlapFile in glob.glob(os.path.join(gregorOutputFolder, name, "index.snp.LD.on.chr*")):
                # if sum(1 for line in open(overlapFile)) > 1:
                df = pandas.read_table(overlapFile, sep='\t', dtype={'inBedPos' : int, 'start' : int, 'end' : int}) 
                mdf = mdf.append(df[['index_SNP','inBedPos','start','end']], ignore_index=True)
                mdf['filename'] = name
        return mdf[['index_SNP','inBedPos','start','end','filename']]
        
    
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='From enrichment in features GREGOR output, fetch the index and/or LD SNP that overlaps a given bed file. Working directory to be the one which contains bed_path_file if it is provided in an argument', usage='python gregorCompileAllReportedOverlappings.py output_T2D/ output.txt -p bedfile.txt')
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
            print("Please provide at least one argument -p or -f")


for files in bedList:
        df = getLdSnps(files, gregorOutputFolder)
        if not df.empty:
                print(files)
                outdf = outdf.append(df, ignore_index=True)

outdf.to_csv(outputfile, sep='\t', index=False, header=['index_SNP','inBedPos','bedStart','bedEnd','filename'])
