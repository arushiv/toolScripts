import pandas
import argparse
import pybedtools
import glob
import os

def getOpts():
    parser = argparse.ArgumentParser(description = 'Format input file for fGWAS', usage = 'python fGWAS_formatInput.py <gwasData> <annotationDirectory> <outputfile> -f chr:pos')
    parser.add_argument('associationFile', help = """input file for association statistics.""")
    parser.add_argument('outputfile', help = """Output filename.""")
    parser.add_argument('-f', '--inputFormat', nargs = '+', help = """format of chrom, pos in input file, followed by column name. Options: chr1:10 colname, 1:10 colname """)
    parser.add_argument('-adir','--annotationDirectory', help = """Directory with annotation bed files.""")
    parser.add_argument('-sa', '--subsetAnnotations', type = str, default="*.bed", help = """Provide wildcard phrase to glob a subset of annotation files from the directory provided. Eg. Islets.*.bed. or *.bed.gz""")
    parser.add_argument('-k', '--keepCols', nargs = '+', help = """Provide header names in current association input that correspond to these in the output: SNPID, Z, F, N""")
    parser.add_argument('-of','--onlyFormat', action='store_const', const="onlyFormat", help="""Only format association file. If the input file is too large, it may be faster to intersect with annotations using bedtools directly rather than within pybedtools. Use this flag in that case.""")
    return parser

def subsetFiles(annotationDirectory, subsetAnnotations):
    outlist = glob.glob(os.path.join(annotationDirectory, subsetAnnotations))
    return outlist

def formatFile(associationFile, keepCols, inputFormat):
    if inputFormat[1] not in keepCols:
        usecols = inputFormat[1] + keepCols
    else:
        usecols = keepCols
    d = pandas.read_csv(associationFile, sep='\t', usecols=usecols)
    if inputFormat[0] == "chr1:10":
        d.loc[:,'CHR'], d.loc[:,'POS'] = d[inputFormat[1]].str.split(":",2).apply(lambda x: x[:2]).str
        d.loc[:,'POS'] = d['POS'].astype(int)
    elif inputFormat[0] == "1:10":
        d.loc[:,'CHR'], d.loc[:,'POS'] = d[inputFormat[1]].str.split(":",2).apply(lambda x: x[:2]).str
        d.loc[:,'POS'] = d['POS'].astype(int)
        d.loc[:,'CHR'] = d['CHR'].map(lambda x: "chr{x}".format(x=x))

    """ Remove weird entries """
    d = d[d['POS'] > 0]

    """ Rename and only retain needed columns """
    # d.loc[:,'start'] = d['POS'] - 1
    neededHeaders = ['SNPID','Z','F','N']
    renameList = dict(zip(keepCols, neededHeaders))
    d.rename(columns = renameList, inplace=True)
    headerlist = ['CHR','POS'] + neededHeaders
    d = d[headerlist]

    """Remove duplicates, retain highest abs(Z)"""
    d.loc[:,'absZ'] = d['Z'].abs()
    d.sort_values(by=['CHR','POS','absZ'], inplace=True)
    d.drop('absZ', axis=1, inplace=True)
    d.drop_duplicates(subset=['CHR','POS'], inplace=True)

    print("assocFileFormatted")
    return d

def getIntersections(assoc_df, fileList):
    assoc_bed = pybedtools.BedTool.from_dataframe(assoc_df).sort()

    for filename in fileList:
        header = os.path.basename(filename).replace(".bed","")
        annot_bed = pybedtools.BedTool(filename)
        assoc_df.loc[:,header] = assoc_bed.intersect(annot_bed, c=True, sorted=True).to_dataframe().iloc[:,-1]
        print(header)
    return assoc_df

if __name__ == '__main__':
    parser = getOpts()
    args = parser.parse_args()
    
    associationFile = args.associationFile
    annotationDirectory = args.annotationDirectory
    outputfile = args.outputfile
    inputFormat = args.inputFormat
    subsetAnnotations = args.subsetAnnotations
    keepCols = args.keepCols
    
    assoc_df = formatFile(associationFile, keepCols, inputFormat)
    
    if args.onlyFormat == "onlyFormat":
        assoc_df.to_csv(outputfile, sep=' ', index=False)
    else:
        fileList = subsetFiles(annotationDirectory, subsetAnnotations)
        dannot = getIntersections(assoc_df, fileList)
        dannot.to_csv(outputfile, sep=' ', index=False)
