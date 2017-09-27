import pandas
import argparse
import pybedtools
import glob
import os

def getOpts():
    parser = argparse.ArgumentParser(description = 'Format input file for fGWAS', usage = 'python fGWAS_formatInput.py <gwasData> <annotationDirectory> <outputfile> -f chr:pos')
    parser.add_argument('associationFile', help = """input file for association statistics.""")
    parser.add_argument('annotationDirectory', help = """Directory with annotation bed files.""")
    parser.add_argument('outputfile', help = """Output filename.""")
    parser.add_argument('-f', '--inputFormat', nargs = '+', help = """format of chrom, pos in input file, followed by column name. Options: chr1:10 colname, 1:10 colname """)

    parser.add_argument('-sa', '--subsetAnnotations', type = str, required = True, help = """Provide wildcard phrase to glob a subset of annotation files from the directory provided. Eg. Islets.*.bed. or *.bed.gz""")
    parser.add_argument('-k', '--keepCols', nargs = '+', help = """Provide header names in current association input that correspond to these in the output: SNPID, Z, F, N""")
    return parser

def subsetFiles(annotationDirectory, subsetAnnotations):
    outlist = glob.glob(os.path.join(annotationDirectory, subsetAnnotations))
    return outlist

def formatFile(associationFile, keepCols, inputFormat):
    d = pandas.read_csv(associationFile, sep='\t')
    if inputFormat[0] == "chr1:10":
        d.loc[:,'CHR'], d.loc[:,'POS'] = d[inputFormat[1]].str.split(":",2).apply(lambda x: x[:2]).str
        d.loc[:,'POS'] = d['POS'].astype(int)
    elif inputFormat[0] == "1:10":
        d.loc[:,'CHR'], d.loc[:,'POS'] = d[inputFormat[1]].str.split(":",2).apply(lambda x: x[:2]).str
        d.loc[:,'POS'] = d['POS'].astype(int)
        d.loc[:,'CHR'] = d['CHR'].map(lambda x: "chr{x}".format(x=x))

    d.loc[:,'start'] = d['POS'] - 1
    neededHeaders = ['SNPID','Z','F','N']
    renameList = dict(zip(keepCols, neededHeaders))
    d.rename(columns = renameList, inplace=True)
    
    headerlist = ['CHR','start','POS'] + neededHeaders
    d = d[headerlist]
    return d

def getIntersections(assoc_df, fileList):
    assoc_bed = pybedtools.BedTool.from_dataframe(assoc_df)

    for filename in fileList:
        header = os.path.basename(filename).replace(".bed","")
        annot_bed = pybedtools.BedTool(filename)
        assoc_df.loc[:,header] = assoc_bed.intersect(annot_bed, c=True).to_dataframe().iloc[:,-1]
    print(assoc_df.head())
        
if __name__ == '__main__':
    parser = getOpts()
    args = parser.parse_args()
    
    associationFile = args.associationFile
    annotationDirectory = args.annotationDirectory
    outputfile = args.outputfile
    inputFormat = args.inputFormat
    subsetAnnotations = args.subsetAnnotations
    keepCols = args.keepCols

    fileList = subsetFiles(annotationDirectory, subsetAnnotations)
    
    assoc_df = formatFile(associationFile, keepCols, inputFormat)    
    dannot = getIntersections(assoc_df, fileList)
    
    # dannot.to_csv(outputfile, sep='\t', index=False)
