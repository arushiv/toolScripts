import pandas
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Merge 2 dataframes by column names', usage='python merge2FilesNyColnames.py file1.dat file2.dat outputfile.dat -c1 col1 col2 col3 -c2 head1 head2 head3')
    parser.add_argument('file1', help="""File 1, tab separated with header names""")
    parser.add_argument('file2', help="""File 2, tab separated with header names""")
    parser.add_argument('outputfile', help ="""Output file.""")
    parser.add_argument('-c1', '--file1colnames', nargs='+', help = """Column headers to merge on for file1""")
    parser.add_argument('-c2', '--file2colnames', nargs='+', help = """Column headers to merge on for file2""")
    parser.add_argument('-t', '--mergetype', default="inner", help = """left, right, outer or inner join. Default='inner'""")

    args = parser.parse_args()
    
    d1 = pandas.read_csv(args.file1, sep='\t')
    d2 = pandas.read_csv(args.file2, sep='\t')
    
    col1 = args.file1colnames
    col2 = args.file2colnames
    mergetype=args.mergetype

    d = pandas.merge(d1,d2, how=mergetype, left_on=col1, right_on=col2) 
    
    d.to_csv(args.outputfile, sep='\t', index=False)
