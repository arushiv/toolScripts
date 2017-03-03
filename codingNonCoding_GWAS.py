# Script function:
# 1. From the given GWAS catalog file, select variants for given traits
# 2. Assemble bed file with uniqe GWAS variants showing associated traits in one column
# 3. Add column indicating if variant lies in coding exon or not using a supplied gencode annotation file


#!/usr/bin/env python
import argparse
import csv
import glob
import os
import pandas
import pybedtools

def gwasSelect(gwasfile, selectTrait):
    dselect = gwasfile.loc[gwasfile['PHENO'].isin(selectTrait)][['CHR','POS','SNP','PHENO']]
    dselect['CHR'] = "chr" + dselect['CHR'].astype(str)
    dunique = dselect.groupby(['CHR','POS','SNP']).agg(lambda x: ",".join(x.tolist())).reset_index()
    dunique['start'] = dunique['POS']
    return pybedtools.BedTool.from_dataframe(dunique[['CHR','start','POS','SNP','PHENO']])


def makeCodingNonCodingGWASBedFile(selectGwasBedFile, gencodefile):
    df = selectGwasBedFile.intersect(gencodefile, c=True).to_dataframe()
    df.rename(columns={'score': 'trait', 'strand': 'overlappingCodingExons'}, inplace=True)
    return df

def printOutput(gwasCodingNonCodingFrame, outputfile):
    gwasCodingNonCodingFrame.to_csv(outputfile, sep='\t', index=False)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="""1. From the given GWAS catalog file, select variants for given traits; add column indicating if variant lies in coding exon or not using a supplied gencode annotation file""")
    parser.add_argument('gwasfile', type=str, help="""gwas catalog file, containing header names as 'CHR', 'PHENO', 'SNP','POS' for columns indicating chrom,phenotype,rsID,position""")
    parser.add_argument('gencodefile', type=str, help="""gencode annotation bed file""")
    parser.add_argument('outputfile', type=str, help="""Output file name. Will contain columns: Chr;Start;end;name;trait;overlappingCodingExons""")
    args = parser.parse_args()

    
    gwasfile = pandas.read_csv(args.gwasfile, sep='\t')
    gencodefile = pybedtools.BedTool(args.gencodefile)
    selectTrait = ["FGlu","FIns","T2D","2hGlu","FGluBMIadj","FInsBMIadj","2hGluBMIadj","Fproinsulin"]
    outputfile = args.outputfile
    
    selectGwasBedFile = gwasSelect(gwasfile, selectTrait)
    gwasCodingNonCodingFrame = makeCodingNonCodingGWASBedFile(selectGwasBedFile, gencodefile)
    printOutput(gwasCodingNonCodingFrame, outputfile)
    
    
