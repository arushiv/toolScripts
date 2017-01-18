# Script function:
# 1. For a given matrix of say, RPKM values for different GTEx samples, calculate mean for samples of each tissue. Samples for each tissue specified in another file
# 2. Using the mean RPKM matrix, compute ESI for each tissue

#!/usr/bin/env python
import argparse
import pandas
import math


def calculateMean(datafile, samplefile):
        sampledf = samplefile.groupby('tissue')['sample'].apply(list)
        sampledf = sampledf[sampledf.apply(lambda x: len(x)) > 25]    # Take those tissues that have more than 25 samples
        
        inputdf = pandas.concat(sampledf.apply(lambda x: chunk.loc[:,x].mean(axis=1)).transpose() for chunk in datafile)

        return inputdf

def log_with_nan(x, y):
    try:
        return math.log(x, y)
    except ValueError:         ## Places with log(0)
        return float('nan')

def entropy_function(x):
    y = -x * log_with_nan(x, 2)
    return y

def filterProteinCoding(d, dCoding):
        filterList = dCoding[dCoding.iloc[:,1]=="protein_coding"].iloc[:,0].tolist()
        return d.ix[filterList]

def calculate_ESI(d):
        d_relative_rpkm = d.div(d.sum(axis=1), axis=0)  ## Relative RPKM = x(g,t)/sum(x(g,t))
        d.loc[:,'entropy'] = d_relative_rpkm.applymap(lambda x: entropy_function(x)).sum(axis=1)    ## Entropy = - sum(relativeRPKM(g,t) * log2(relative RPKM(g,t)))
        d_qvalues = d_relative_rpkm.applymap(lambda x: log_with_nan(x, 2)).apply(lambda x: d.loc[:,'entropy'] - x, axis=0)   ## Q(g,t) = Entropy(g) - log2(relativeRPKM(g,t))
        d_max_q = d_qvalues.max(axis=0).to_frame().transpose()        ## maxQ = max(Q(t))

        d_esi = d_qvalues.apply(lambda x: 1 - x/d_max_q.squeeze(), axis=1)
        return d_esi

    
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='For a given matrix of say, RPKM values for different GTEx samples, calculate mean RPKM for each tissue. Sample IDs corresponding each tissue are specified in another file. Then calculate ESI for each tissue', usage='python esiScoreAfterMean.py datafile.txt samplefile.txt output.txt')
    parser.add_argument('datafile', help="""Input Matrix. Tab separated, with header specifying sample ID. First two columns taken as 'gene name' and 'description'. Rest columns are to be sample IDs. Designed to handle large GTEx data frames, this file is read in chunks of 5000 lines""")
    parser.add_argument('samplefile', help="""Tab separated file with 2 columns with headers: 'tissue' 'sample', specifying sample IDs corresponding to each tissue type.""")
    parser.add_argument('outputfile', help ="""Output file.""")
    parser.add_argument('-f', '--geneFilterList', default = 'NoFileSupplied', help = """File with gene ID as first column and gene type as second column. No header. will filter for genes with type = "protein_coding" before calculating ESI""")
    args = parser.parse_args()

    datafile = pandas.read_csv(args.datafile, sep='\t', index_col=[0,1], comment='#', header=1, low_memory=False, chunksize = 5000, iterator=True)
    samplefile = pandas.read_csv(args.samplefile, sep='\t')
    outputfile = args.outputfile
    # Calculate Mean for each tissue given the sample names
    inputdf = calculateMean(datafile, samplefile)

    if args.geneFilterList != "NoFileSupplied":
            dCoding = pandas.read_csv(args.geneFilterList, sep='\t', header=None)
            inputdf = filterProteinCoding(inputdf, dCoding)
        
    # Calculate ESI
    d_esi = calculate_ESI(inputdf)
    d_esi.to_csv(outputfile, sep='\t', na_rep="NA")
