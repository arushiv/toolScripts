import pandas
import argparse
import pybedtools
import glob
import os
import functools

def getOpts():
    parser = argparse.ArgumentParser(description = 'Merge multiple files with some columns of one main file', usage = 'python mergeMultipleFiles.py <mainFile> <> <outputfile> -f chr:pos')
    parser.add_argument('--mainFile', required = True, help = """Main input files with columns that all files have in common.""")
    parser.add_argument('--outputfile', required = True, help = """Output filename.""")
    parser.add_argument('-f', '--filesToMerge', nargs = '+', help = """Space separated paths of files to merge""")
    parser.add_argument('-c1', '--colnamesMainFile', nargs = '+', help = """Provide header names if mainFile doesn't have a header""")
    parser.add_argument('-c2', '--colnamesMergeFile', nargs = '+', help = """Provide header names if filesToMerge don't have headers""")
    parser.add_argument('-on','--mergeOn', nargs = '+', required = True, help = """Column names to merge on.""")
    parser.add_argument('--renameFilename', help = """Rename the last column of filesToMerge by own filename. Used in eg. input formatting for fGWAS.""")

    return parser

def readNReplace(filename, colnamesMergeFile, renameFilename):
    if colnamesMergeFile is not None:
        df = pandas.read_csv(filename, sep='\t', header=None, names=colnamesMergeFile)
    else:
        df = pandas.read_csv(filename, sep='\t')
        
    headerName = os.path.basename(filename).replace(".bed","")
    df.rename(columns={renameFilename : headerName}, inplace=True)
    return df
        
def read(filename, colnamesMergeFile):
    if colnamesMergeFile is not None:
        df = pandas.read_csv(filename, sep='\t', header=None, names=colnamesMergeFile)
    else:
        df = pandas.read_csv(filename, sep='\t')
    return df


if __name__ == '__main__':
    parser = getOpts()
    args = parser.parse_args()

    outputfile = args.outputfile
    colnamesMainFile = args.colnamesMainFile
    colnamesMergeFile = args.colnamesMergeFile
    renameFilename = args.renameFilename
    
    if colnamesMainFile is not None:
        dmain = pandas.read_csv(args.mainFile, sep='\t', header=None, names=args.colnamesMainFile)
    else:
        dmain = pandas.read_csv(args.mainFile, sep='\t')
    print("mainFile read")
    
    if renameFilename is not None:
        dflist = [readNReplace(filename, colnamesMergeFile, renameFilename) for filename in args.filesToMerge]
    else:
        dflist = [read(filename, colnamesMergeFile) for filename in args.filesToMerge]
    print(dflist)
    
    df_final = functools.reduce(lambda left,right: pandas.merge(left, right, on=args.mergeOn), dflist)
    print(df_final.head())
    df_final.to_csv(args.outputfile, sep='\t', index=False, compression="gzip")
