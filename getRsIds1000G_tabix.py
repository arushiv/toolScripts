import pandas
import argparse
import glob
import subprocess as sp
import os



def tabixCommand(row, vcfFileDirectory):
    cmd = "tabix -fh %s %s | grep -v '#' | cut -f3 | grep -v esv | sort | uniq" %(os.path.join(vcfFileDirectory, "ALL.chr%s.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz" %(row['SNPchr'])), row['snp1'])
    return sp.check_output(cmd, shell=True).replace("\n","")

def useTabix(df, vcfFileDirectory):
    df.loc[:,'snpRsID'] = df.apply(lambda row: tabixCommand(row, vcfFileDirectory), axis=1)
    df.drop(['snp1','SNPchr','SNPEnd'], axis=1, inplace=True)
    return df


    
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='For a bed file with chrom and pos fields for SNPs, fetch rsids and output in a dataframe', usage='python getRsIds1000G.py inputfile -o outputfile')
    parser.add_argument('inputfile', type=str, help="""The input file. Should contain one colum with header 'snp' and values formatted as 'chr1:10000' """)
    parser.add_argument('-v', '--vcfFileDirectory', type=str, default = "/lab/data/genomes/human/hg19/1000GenomesDownloads/", help="""The directory with vcf files. Default = `/lab/data/genomes/human/hg19/1000GenomesDownloads/`""")
    parser.add_argument('outputfile', help="""output file name""")

    args = parser.parse_args()
    vcfFileDirectory = args.vcfFileDirectory
    df = pandas.read_csv(args.inputfile, sep='\t')

    df.loc[:,'SNPchr'], df.loc[:,'SNPEnd'] = df['snp'].str.split(':').str
    df.loc[:,'SNPchr'] = df['SNPchr'].str.replace("chr","").astype(int)

    def organize(x):
        return "%s:%s-%s" %(x['SNPchr'], x['SNPEnd'], x['SNPEnd'])
    
    df.loc[:,'snp1'] = df.apply(lambda x: organize(x), axis=1) 
           
    ndf = useTabix(df, vcfFileDirectory)
    # ndf = pandas.concat([readVcfAndMerge(df, vcf) for vcf in filelist])
    
    ndf.to_csv(args.outputfile, index=False, sep='\t', na_rep="NA")
