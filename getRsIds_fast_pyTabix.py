import pandas
import tabix
import numpy
import argparse
import time
import sys

def getOpts():
    parser = argparse.ArgumentParser(description='For an input file with chrom and pos fields for SNPs, fetch rsids and output in a dataframe. '
                                     'Uses dbsnp 150 vcf by default '
                                     '!!!IMPORTANT!!!: The script will return the first rsID that overlaps the SNP position, this might not be what is needed if there are overlapping indels. '
                                     'If you only need SNPs, subset your vcf first using vcftools --remove-indels and then use this script.' 
                                     'This script uses chunks and is fast for, say n hundred thousand queries, but could take hours for example working with GWAS summary data with many millions of quesries.',
                                     usage='python getRsIds_fast_pyTabix.py --inputfile <inputfile> --outputfile <outputfile> --vcfFile <filepath> --chrom chromosome --pos position')
    parser.add_argument('--inputfile', required=True, type=str,
                        help="""Tab delimited input file. Should contain at least chrom and pos information. See other options to specify format""")
    parser.add_argument('-v', '--vcf', type=str, default = "/lab/data/reference/human/hg19/annot/dbsnp150_variants/common_all_20170710.vcf.gz",
                        help="""The vcf file. Default = `/lab/data/reference/human/hg19/annot/dbsnp150_variants/common_all_20170710.vcf.gz.`""")
    parser.add_argument('--header-pos', type=int,
                        help="""Row number (0 indexed) which is to be taken as header. If no header, specify the --header-names argement instead.""")
    parser.add_argument('--header-names', type=str, nargs='+',
                        help="""If input file in a headerless, provide column names. """)
    parser.add_argument('-c', '--chrom', required=True, type=str,
                        help="""Column name for chromosome """)
    parser.add_argument('-p', '--pos', required=True, type=str,
                        help="""Column name for SNP position """)
    parser.add_argument('--chunksize', type=int, default=300000,
                        help="""Number of rows in a chunk as inputfile is read in chunks """)
    parser.add_argument('--chrom-type', 
                        help="""type of elements in the chrom column. str if column contains 1, 3, X, Y; int if column contains only chrom numbers""")
    parser.add_argument('--subset-cols', type=str, nargs='+',
                        help="""If input file is too large, subsetting to read in only required columns can be helpful. 
                        Provide these column names.""")

    parser.add_argument('--outputfile', required=True,
                        help="""output file name. Will contain an extra column named 'snpRsID'.""")
    args = parser.parse_args()
    return args

def peek(iterable):
    """Check if tabix output has some values. Return the third value that is the rsID. If the tabix output is empty, return nan """
    try:
        value = next(iterable)[2]
    except StopIteration:
        return numpy.nan
    except TypeError:
        return numpy.nan
    return value

def query_function(chrom, pos, chromType, opened_tabix):
    """If some weird chromosome is not found in the vcf, tabix throws TabixError. Just return nan in
    such cases so rest of the positions can be queried"""
    function_dictionary = {
        'str' : opened_tabix.query,
        'int' : opened_tabix.queryi,
    }
    try:
        value = function_dictionary[chromType](chrom, pos, pos)
    except tabix.TabixError:
        return numpy.nan
    return value
    
def apply_on_chunk(chunk, chromType, opened_tabix):
    chunk.loc[:,'snpRsID'] = chunk.apply(lambda x: peek(query_function(x['chrom'], x['pos'], chromType, opened_tabix)), axis=1)
    print("chunk processed")
    return chunk
  
    

if __name__ == '__main__':

    args = getOpts()
    if ((args.header_pos is not None and args.header_names is not None) or ( args.header_pos is None and args.header_names is None)):
        print("Provide either --header-pos or --header-names argument")
        sys.exit()

    d = pandas.read_csv(args.inputfile, sep='\t', header=args.header_pos, names=args.header_names, usecols=args.subset_cols, chunksize=args.chunksize)
    
    tb = tabix.open(args.vcf)

    l = [apply_on_chunk(chunk, args.chrom_type, tb) for chunk in d]
   
    pandas.concat(l).to_csv(args.outputfile, index=False, sep='\t', na_rep="NA")


