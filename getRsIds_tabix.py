import pandas
import argparse
import subprocess as sp


def useTabix(df, vcfFile):
    def tabixCommand(row, vcfFile):
        cmd = r"""tabix {vcfFile} {region} | cut -f3 | grep -v esv | sort | uniq""".format(vcfFile=vcfFile, region=row['snpPosToFetch'])
        return sp.check_output(cmd, shell=True)

    df.loc[:,'snpRsID'] = df.apply(lambda row: tabixCommand(row, vcfFile), axis=1).str.decode('ascii').str.strip()
    df.drop(['snpPosToFetch'], axis=1, inplace=True)
    return df

def fixInputFile(filename, header, chrom_col, pos_col):
    if header is not None:
        df = pandas.read_csv(filename, sep='\t', header=None, names=header)
    else:
        df = pandas.read_csv(filename, sep='\t')

    def organize(x):
        out = "{chrom}:{pos}-{pos}".format(chrom=x[chrom_col], pos=x[pos_col])
        return out
    
    df.loc[:,'snpPosToFetch'] = df.apply(lambda x: organize(x), axis=1) 
    df.loc[:,'snpPosToFetch'] = df['snpPosToFetch'].str.replace("chr","")

    return df
    
def getOpts():
    parser = argparse.ArgumentParser(description='For a input file with chrom and pos fields for SNPs, fetch rsids and output in a dataframe',
                                     usage='python getRsIds_tabix.py --inputfile <inputfile> --outputfile <outputfile> --vcfFile <filepath> --chrom chromosome --pos position')
    parser.add_argument('--inputfile', required=True, type=str,
                        help="""Tab delimited input file. Should contain at least chrom and pos information. See other options to specify format""")
    parser.add_argument('-v', '--vcfFile', type=str, default = "/lab/data/reference/human/hg19/annot/dbsnp150_variants/common_all_20170710.vcf.gz",
                        help="""The vcf file. Default = `/lab/data/reference/human/hg19/annot/dbsnp150_variants/common_all_20170710.vcf.gz`""")
    parser.add_argument('--header', type=str, nargs='+',
                        help="""If input file in a headerless BED format, provide column names. """)
    parser.add_argument('-c', '--chrom', required=True, type=str,
                        help="""Column name for chromosome """)
    parser.add_argument('-p', '--pos', required=True, type=str,
                        help="""Column name for SNP position """)
    parser.add_argument('--outputfile', required=True,
                        help="""output file name. Will contain an extra column named 'snpRsID'.""")
    args = parser.parse_args()
    return args
    
if __name__ == '__main__':

    args = getOpts()

    df = fixInputFile(args.inputfile, args.header, args.chrom, args.pos)
    print(df)
    ndf = useTabix(df, args.vcfFile)
    
    ndf.to_csv(args.outputfile, index=False, sep='\t', na_rep="NA")
