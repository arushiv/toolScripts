import pandas
import argparse


def getOpts():
    parser = argparse.ArgumentParser(description='Script to split file by factors in a column. Output files will be saved as <givePhrase>.<factorName>.<bed>', usage='python splitFileByColumnFactors.py <inputfile> -c colname -o outputfilePhrase ')
    parser.add_argument('inputfile', help="""input file, tab separated""")
    parser.add_argument('-c', '--colname', nargs='+', help="""name of column header""")
    parser.add_argument('-ih', '--inputFileHeader', nargs='+', help = """If the input file doesn't have headers, provide column names""")
    parser.add_argument('-o', '--outputFilePhrase', help="""phrase of output file""")
    parser.add_argument('-k', '--keepCols', nargs='+', help = """Keep these columns in the final output.""")


    return parser

def splitfile(inputfile, colname, outputFilePhrase, inputFileHeader, keepCols):
    if inputFileHeader is not None:
        d = pandas.read_csv(inputfile, sep='\t', header=None, names=inputFileHeader)
    else:
        d = pandas.read_csv(inputfile, sep='\t')
    for name, group in d.groupby(colname):
        filename = "{outputFilePhrase}.{name}.bed".format(outputFilePhrase=outputFilePhrase, name=name)
        if keepCols is not None:
            group[keepCols].to_csv(filename, sep='\t', index=False)
        else:
            group.to_csv(filename, sep='\t', index=False)
    

if __name__ == '__main__':
    parser = getOpts()
    args = parser.parse_args()
    inputfile = args.inputfile
    colname = args.colname
    outputFilePhrase = args.outputFilePhrase

    splitfile(inputfile, colname, outputFilePhrase, args.inputFileHeader, args.keepCols)

