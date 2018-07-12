import pandas
import tabix
import numpy
import argparse

def getOpts():
    parser = argparse.ArgumentParser(description='For an input file with chrom and pos fields for SNPs, fetch rsids and output in a dataframe. '
                                     'Uses dbsnp 150 vcf by default '
                                     '!!!IMPORTANT!!!: The script will return the first rsID that overlaps the SNP position, this might not be what is needed if there are overlapping indels. '
                                     'If you only need SNPs, subset your vcf first using vcftools --remove-indels and then use this script. ',
                                     usage='python getRsIds_fast_pyTabix.py --inputfile <inputfile> --outputfile <outputfile> --vcfFile <filepath> --chrom chromosome --pos position')
    parser.add_argument('--inputfile', required=True, type=str,
                        help="""Tab delimited input file. Should contain at least chrom and pos information. See other options to specify format""")
    parser.add_argument('-v', '--vcfFile', type=str, default = "/lab/data/reference/human/hg19/annot/dbsnp150_variants/common_all_20170710.vcf.gz",
                        help="""The vcf file. Default = `/lab/data/reference/human/hg19/annot/dbsnp150_variants/common_all_20170710.vcf.gz.`""")
    parser.add_argument('--header', type=str, nargs='+',
                        help="""If input file in a headerless BED format, provide column names. """)
    parser.add_argument('-c', '--chrom', required=True, type=str,
                        help="""Column name for chromosome """)
    parser.add_argument('-p', '--pos', required=True, type=str,
                        help="""Column name for SNP position """)
    parser.add_argument('--chrType', required=True, type=str,
                        help="""write either "str" or "int" to specify if the chrom column has values of the format chr10 or 10 """)
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
    return value

if __name__ == '__main__':

    args = getOpts()
    
    if args.header is not None:
        d = pandas.read_csv(args.inputfile, sep='\t', header=None, names=args.header)
    else:
        d = pandas.read_csv(args.inputfile, sep='\t')

    tb = tabix.open(snakemake.input.vcf)

    """ Get rsid """
    if args.chromType == "str":
        d.loc[:,'snpRsID'] = d.apply(lambda x: peek(tb.query(x[args.chrom], x[args.pos], x[args.pos])), axis=1)

    elif args.chromType == "int":
        d.loc[:,'snpRsID'] = d.apply(lambda x: peek(tb.queryi(x[args.chrom], x[args.pos], x[args.pos])), axis=1)

    d.to_csv(args.outputfile, index=False, sep='\t', na_rep="NA")


