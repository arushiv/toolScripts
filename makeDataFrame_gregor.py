# Run GREGOR: Fix conf files, submit drmr commands, compile dataframe and plot 
#!/usr/bin/env python


import argparse
import pandas
import glob
                   
def makedataframe(filename, fieldSeparator, outputfilename):
    job_info = filename.replace('/StatisticSummaryFile.txt','').replace('output_','')

    df = pandas.read_csv(filename, sep='\t')
    df['Bed_File'] = df['Bed_File'].str.replace('.bed','')
    df['job_info'] = job_info
    splitBeddf = df['Bed_File'].str.split(fieldSeparator, expand=True)
    splitJobdf = df['job_info'].str.split("_", expand=True)
    splitSnpdf = splitJobdf.iloc[:,0].str.split(".", expand=True)
    outdf = pandas.concat([splitSnpdf, splitJobdf.drop(splitJobdf.columns[0], 1), splitBeddf, df.drop(['Bed_File','job_info'], 1)], 1)
    
    outdf.to_csv(outputfilename, mode='a', header=False, index=False, sep='\t', na_rep="NA")
    
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Compile dataframe from GREGOR output folders.')
    parser.add_argument('filename', help="""The StatisticSummaryFile from GREGOR.""")
    parser.add_argument('-f','--nameFieldSeparator', type=str, default='.', help="""Field separator to make columns from bed file name. (Default='.')""")
    parser.add_argument('outputfilename', help="""Output file name.""")
    args = parser.parse_args()

    filename = args.filename
    fieldSeparator = args.nameFieldSeparator
    outputfilename = args.outputfilename

    makedataframe(filename, fieldSeparator, outputfilename)
