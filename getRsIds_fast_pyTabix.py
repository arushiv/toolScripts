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
                        help="""type of chrom column. object or int64 """)
    parser.add_argument('--subset-cols', type=str, nargs='+',
                        help="""If input file is too large, subsetting to read in only required columns can be helpful. 
                        Provide these column names.""")

    parser.add_argument('--outputfile', required=True,
                        help="""output file name. Will contain an extra column named 'snpRsID'.""")
    args = parser.parse_args()
    return args

def myQuery_int(chrom, pos):
    """If some weird chromosome is not found in the vcf, tabix throws TabixError. Just return nan in
    such cases so rest of the positions can be queried"""
    try:
        value = tb.queryi(chrom, pos, pos)
    except tabix.TabixError:
        return numpy.nan
    return value

def myQuery_str(chrom, pos):
    """
    If chrom column is of the datatype object, i.e. has str values:
    If some weird chromosome is not found in the vcf, tabix throws TabixError. Just return nan in
    such cases so rest of the positions can be queried"""
    try:
        value = tb.query(chrom, pos, pos)
    except tabix.TabixError:
        return numpy.nan
    return value

def peek(iterable):
    """Check if tabix output has some values. Return the third value that is the rsID. If the tabix output is empty, return nan """
    try:
        value = next(iterable)[2]
    except StopIteration:
        return numpy.nan
    except TypeError:
        return numpy.nan
    return value

  
    

if __name__ == '__main__':

    args = getOpts()

    if (not args.header-pos and not args.header-names) or ( args.header-pos and args.header-names):
        print("Provide either --header-pos or --header-names argument")
        sys.exit()

    d = pandas.read_csv(args.inputfile, sep='\t', header=args.header-pos, names=args.header-names, usecols=args.subsetCols, chunksize=args.chunksize)
    
    tb = tabix.open(args.vcf)

    # chromType = d[[args.chrom]].dtypes[0]
    # def apply_on_chunk(chunk):
    #     chunk.loc[:,'snpRsID'] = chunk.apply(lambda x: peek(myQuery_str(x['chrom'], x['pos'])), axis=1)
    #     print("chunk processed")
    #     return chunk
    
    """ Get rsid """
    if args.chromType == "object":
        # l = [apply_on_chunk(chunk) for chunk in d]
        # tempcol = "tempcol{time}".format(time = time.time())
        # d.loc[:, tempcol] = d[args.chrom].str.replace("chr", "")
        for chunk in d:
            chunk.loc[:,'snpRsID'] = chunk.apply(lambda x: peek(myQuery_str(x['chrom'], x['pos'])), axis=1)
            l.append(chunk)
            print("chunk processed")
        
        # # d.loc[:,'snpRsID'] = d.apply(lambda x: peek(myQuery_str(x[tempcol], x[args.pos])), axis=1)
        # d.drop(tempcol, axis=1, inplace=True)
        
    elif args.chromType == "int64":
        for chunk in d:
            chunk.loc[:,'snpRsID'] = chunk.apply(lambda x: peek(myQuery_int(x['chrom'], x['pos'])), axis=1)
            l.append(chunk)

        # d.loc[:,'snpRsID'] = d.apply(lambda x: peek(myQuery_int(x[temp_chrom_column], x[args.pos])), axis=1)


            
    # if args.chrType == "str":
    #     d.loc[:,'snpRsID'] = d.apply(lambda x: peek(tb.query(x[args.chrom], x[args.pos], x[args.pos])), axis=1)

    # elif args.chrType == "int":
    #     d.loc[:,'snpRsID'] = d.apply(lambda x: peek(tb.queryi(x[args.chrom], x[args.pos], x[args.pos])), axis=1)

   
    pandas.concat(l).to_csv(args.outputfile, index=False, sep='\t', na_rep="NA")


