import pandas
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Subest file by column names', usage='python subsetFileNyColName.py file.dat -c col1 col2 col3 outfile.dat')
    parser.add_argument('file1', help="""File 1, with header names""")
    parser.add_argument('outputfile', help ="""Output file.""")
    parser.add_argument('-c', '--filecolnames', nargs='+', help = """Column headers to subset file on""")
    parser.add_argument('-sep', '--separatorfile', default="\t", help = """File 1 field separator. Default = '\\t'""")

    args = parser.parse_args()
    
    d = pandas.read_csv(args.file1, sep=args.separatorfile)
    
    col = args.filecolnames

    dout = d[col]
        
    dout.to_csv(args.outputfile, sep='\t', index=False)
