# Run GREGOR: Fix conf files, submit drmr commands, compile dataframe and plot 
#!/usr/bin/env python
# 
# from __future__ import print_function
import argparse
import collections
import csv
import gzip
import magic
import pprint
import re
import operator
import glob
import math
import sys
import subprocess as sp
import os
import itertools
import shutil
import fileinput
import fnmatch

# class RunGregor(object):
#     def __init__(self):
#         pass
    
def makedataframe(filename, fieldSeparator, outputfilename):
    with open(outputfilename, 'a') as fout:
        with open(filename) as f:
            for line in f:
                if "Bed_File" not in line:
                    s=line.split(' ')
                    snpcolumns = '\t'.join(((filename.split('/')[0]).replace('output_','')).split('_'))
                    namecolumns = '\t'.join((s[0].replace('.bed','')).split(fieldSeparator))
                    datacolumn = '\t'.join(s[1:])
                    fout.write("%s\t%s\t%s"%(snpcolumns,namecolumns,datacolumn))
                   
    
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Run GREGOR: Fix conf files, submit drmr commands, compile dataframe and plot.')
    parser.add_argument('filename', help="""The StatisticSummaryFile from GREGOR.""")
    parser.add_argument('-f','--nameFieldSeparator', type=str, default='.', help="""Field separator to make columns from bed file name. (Default='.')""")
    parser.add_argument('outputfilename', help="""Output file name.""")
    args = parser.parse_args()

    filename=args.filename
    fieldSeparator=args.nameFieldSeparator
    outputfilename=args.outputfilename

    makedataframe(filename, fieldSeparator, outputfilename)
