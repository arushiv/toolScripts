import pandas
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert excel file to dataframe',
                                     usage='python convertExcelToDataframe.py inputfile.xlsx outputfile.dat -c col1 col2 col3')
    parser.add_argument('file1', help="""File 1, with header names""")
    parser.add_argument('outputfile', help ="""Output file.""")
    parser.add_argument('-c', '--filecolnames', nargs='+', help = """Column headers to subset file on""")
    parser.add_argument('-sep', '--separatorfile', default="\t", help = """Output file field separator. Default = '\\t'""")
    parser.add_argument('--sheet', help = """Sheet name to read.""")
    args = parser.parse_args()

    if args.sheet is not None:
        d = pandas.read_excel(args.file1, sheet_name=args.sheet)
    else:
        d = pandas.read_excel(args.file1)
    
    col = args.filecolnames

    if col is not None:
        d = d[col]
    
    d.to_csv(args.outputfile, sep=args.separatorfile, index=False)
