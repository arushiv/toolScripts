# Match rsID on provided column


#!/usr/bin/env python


#from __future__ import print_function
import argparse
import csv
import pandas

# def makedummy1(row):
#     return "0"

# def makedummy2(row):
#     return "."

# def makedummy3(row):
#     return "1"

# def assignColour(row):
#     if row['buddyrs'] == row['bestrs']:
#         return '75,0,130'
#     elif row['buddyrs'] != row['bestrs']:
#         return '252,58,47'

def print_missing_variants(filebuddies, filebesteqtl):
    newtable = pandas.merge(filebuddies, filebesteqtl, how='outer', on='snpcol')
    return newtable


def print_to_file(newtable, outputfile):
                # self.header.append(self.pval)
        newtable.to_csv(outputfile, sep='\t', index=False, na_rep="NaN")

if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('filebuddies', help="""The file with LD buddies for each best eQTL. Originally designed for header 'chr start pos bestrs buddyrs'.""")
    parser.add_argument('filebesteqtl', help="""The file with gene name and information for each best eQTL. Orignially designed for header 'best_chr best_start best_end best_allele geneID beta direction genetype gene bestrs'.""")
    # parser.add_argument('-f', '--column_buddyFile', type=str, default="bestrs", help="""The column number in the LD buddy file to match.""")
    # parser.add_argument('-qf', '--column_bestEqtlFile', type=str, default="bestrs", help="""The column number in best eQTL file.""")
    parser.add_argument('outputfile', help="""Name of the output file. Originally designed to output 'chr start pos buddyrs 0 . start end colour gene direction' for a bedDetail format""")

    args = parser.parse_args()

    filebuddies = pandas.read_csv(args.filebuddies, sep='\t')
    filebesteqtl = pandas.read_csv(args.filebesteqtl, sep='\t')
    # column_buddyFile = args.column_buddyFile
    outputfile = args.outputfile
    # column_bestEqtlFile = args.column_bestEqtlFile
    
    # print_to_file(makeTrackBed(filebuddies, filebesteqtl, column_buddyFile, column_bestEqtlFile), outputfile)
    print_to_file(print_missing_variants(filebuddies, filebesteqtl), outputfile)
    # print_missing_variants(filebuddies, filebesteqtl)
