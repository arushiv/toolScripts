# Standard python template for opening files 
#!/usr/bin/env python


#from __future__ import print_function
import argparse
import collections
import csv
import gzip
import magic
import pprint
import re
import operator
import glob

def open_maybe_gzipped(f):
    with open(f) as test_read:
        mimetype = magic.from_buffer(test_read.read(1024), mime=True).decode('utf-8')
    if mimetype == 'application/x-gzip':
        f = gzip.open(f, mode='rt') # Python3: , encoding='utf-8')
    else:
        f = open(f, 'rt')
    return f

def matchColor(bedfile, colorfile, bedcolumn, outfile, statecolumn, colorcolumn, statenumbercolumn):
    with open_maybe_gzipped(bedfile) as bfile, open_maybe_gzipped(colorfile) as cfile, open(outfile, 'w') as ofile:
        cfile.readline()
        for stuff in cfile:
            l = stuff.split('\t')
            state = "%s_%s"%(l[statenumbercolumn-1],l[statecolumn-1].replace("/","_"))
            color = l[colorcolumn-1].replace("\n","")
            
            for line in bfile:
                lb = line.split('\t')
                lb[bedcolumn-1] = lb[bedcolumn-1].replace("\n","")
                if lb[bedcolumn-1] == state:
                    lb.append(color)
                    ofile.write("%s\n"%("\t".join(lb)))
            bfile.seek(0)    

if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('bedfile', help="""The bed file with nth column as the state name (default n = 4th column)""")
    parser.add_argument('colorfile', help="""The file with colors for each state""")
    parser.add_argument('-f','--statefieldbedfile', nargs='?', type=int, default=4, help="""The column number in the bed file containing the state name (default=4)""")
    parser.add_argument('-cs','--statenamefieldcolorfile', nargs='?', type=int, default=5, help="""The column number in the color file containing the state name (default=5)""")
    parser.add_argument('-csn','--statenumberfieldcolorfile', nargs='?', type=int, default=3, help="""The column number in the color file containing the state number (default=3)""")
    parser.add_argument('-cc','--colorfieldcolorfile', nargs='?', type=int, default=7, help="""The column number in the color file containing the state name (default=7)""")
    parser.add_argument('outputfile', help="""Name of the output file""")

    args = parser.parse_args()

    bedfile = args.bedfile
    colorfile = args.colorfile
    bedcolumn = args.statefieldbedfile
    outfile = args.outputfile
    statecolumn = args.statenamefieldcolorfile
    colorcolumn = args.colorfieldcolorfile
    statenumbercolumn = args.statenumberfieldcolorfile

    matchColor(bedfile, colorfile, bedcolumn, outfile, statecolumn, colorcolumn, statenumbercolumn)
                
