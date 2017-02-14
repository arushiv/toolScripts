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
                        

def makeLdBedFile(df):  # From index SNP and LD buddy dataframe, make bed file with chrom start end for each snp and LD snp

        s = df['LD_buddy_pos'].str.split("|", expand=True).apply(pandas.Series, 1).stack(dropna=True)
        s.index = s.index.droplevel(-1)
        s.name = 'pos'
        dreturn = df.join(s)
        dreturn["chrom"], dreturn["indexSnp"] = zip(*dreturn["index_SNP"].str.split(':').tolist())
        
        del dreturn["index_SNP"]
        dreturn = dreturn[["chrom","pos","pos","indexSnp"]]
        dreturn.iloc[:,1] = dreturn.iloc[:,1].apply(int) - 1

        return dreturn


    
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Compile LD SNPs and index SNPs themselves form the index_SNP folder from a GREGOR output.', usage='python gregorFindAllLdSnps.py <path/to/index_SnpFile/index.snp.LD.txt>')
    parser.add_argument('gregorOutputFile', type=str, help="""GREGOR output file for index.snp.LD.txt """)
    parser.add_argument('outputfile', help ="""Output file. Will contain these new columns: chrom_overlapSNP start_overlapSNP end_overlapSNP pos_indexSNP""")
    args = parser.parse_args()

    outputfile = args.outputfile
    outdf = pandas.DataFrame()
    indexSNPFile = pandas.read_csv(args.gregorOutputFile, sep='\t')
    makeLdBedFile(indexSNPFile).to_csv(outputfile, sep='\t', index=False, header=False)
    # if not df.empty:
#             ldBedFile = makeLdBedFile(df)
#             outdf = outdf.append(ldBedFile, ignore_index=True)

# outdf.to_csv(outputfile, sep='\t', index=False, header=False)
