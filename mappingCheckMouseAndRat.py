import pandas
import numpy
import argparse

def mergeFiles(main, mouseMapped, ratMapped):
    main_mouse = pandas.merge(main, mouseMapped, how="left", on="snpinfo")
    main_mouse_rat = pandas.merge(main_mouse, ratMapped, how="left", on="snpinfo")
    main_mouse_rat = main_mouse_rat[['snpinfo','me','re']]
    
    main_mouse_rat.loc[:,'mapping'] = numpy.where((main_mouse_rat.me.notnull() & main_mouse_rat.re.notnull()), "MouseAndRat", "_")
    main_mouse_rat.loc[:,'mapping'] = numpy.where((main_mouse_rat.me.notnull() & main_mouse_rat.re.isnull()), "onlyMouse", main_mouse_rat['mapping'])
    main_mouse_rat.loc[:,'mapping'] = numpy.where((main_mouse_rat.me.isnull() & main_mouse_rat.re.notnull()), "onlyRat", main_mouse_rat['mapping'])
    main_mouse_rat.loc[:,'mapping'] = numpy.where((main_mouse_rat.me.isnull() & main_mouse_rat.re.isnull()), "NotMappable", main_mouse_rat['mapping'])

    main_mouse_rat.loc[:,'chrom'], main_mouse_rat.loc[:,'snpPos'], main_mouse_rat.loc[:,'FGluIndex'], main_mouse_rat.loc[:,'T2DIndex'] = main_mouse_rat.snpinfo.str.split(":").str

    del main_mouse_rat['snpinfo']
    return main_mouse_rat
# main_mouse_rat.loc[:,'proxy'] = main_mouse_rat[['chrom','proxy']].apply(lambda x: '_'.join(x), axis=1)
# main_mouse_rat.loc[:,'index'] = main_mouse_rat[['chrom','index']].apply(lambda x: '_'.join(x), axis=1)
# del main_mouse_rat['chrom']


def saveFile(df, outputfile):
    df.to_csv(outputfile,  sep='\t', na_rep="NA", index=False)



if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='From an hg19 bedfile (to bnMap) and bnMapper outputs for mapping to mouse and rat, make one dataframe noting SNPs that map to mouse and/or rat.', usage='python mappingCheckMouseAndRat.py inputfile mousefile ratfile outputfile')
    parser.add_argument('inputfile', type=str, help="""SNP bed file with chrom, start, end, snpinfo. No header""")
    parser.add_argument('mousefile', type=str, help="""SNP bed file with chrom, start, end, snpinfo. No header""")
    parser.add_argument('ratfile', type=str, help="""SNP bed file with chrom, start, end, snpinfo. No header""")
    parser.add_argument('outputfile', type=str, help="""Outputfile name""")
    parser.add_argument('-s', '--snpinfosplit', nargs='+', default=['chrom','snpPos','FGluIndex','T2DIndex'], help="""Supply comma separated list of names into which the snpinfo column should be split into. Default = ['chrom','snpPos','FGluIndex','T2DIndex']""") 

    
    args = parser.parse_args()


    main = pandas.read_csv(args.inputfile, sep='\t', header=None, names = ['c','s','e','snpinfo'])
    mouseMapped = pandas.read_csv(args.mousefile, sep='\t', header=None, names = ['mc','ms','me','snpinfo'])
    ratMapped = pandas.read_csv(args.ratMapped, sep='\t', header=None, names = ['rc','rs','re','snpinfo'])


    df = mergeFiles(main, mouseMapped, ratMapped)

    saveFile(df, args.outputfile)

