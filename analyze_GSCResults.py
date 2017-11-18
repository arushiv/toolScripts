#!/usr/bin/env python
import argparse
import pandas
import math
import glob
import os

def log_with_nan(x, y):
    try:
        return math.log(x, y)
    except ValueError:         ## Places with log(0)
        return float('nan')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Make dataframe from GSC results.', usage='python ~arushiv/toolScripts/analyze_GSCResults.py 100 -s .dat -d outputDir/ ')
    parser.add_argument('-s','--inputString', type=str, default='*.txt', help="""Files which should be parsed. (default: *.txt)""")
    parser.add_argument('-d','--resultDir', type=str, default='.', help="""Directory where result files reside. (default: current Directory)""")
    parser.add_argument('runs', type=int, help="""Provide how many runs was GSC run for.""")
    parser.add_argument('outputfile', type=str, help="""Output file name.""")
    args = parser.parse_args()

    path = r'.'                     # use your path
    resultDir = args.resultDir
    runs = args.runs
    inputString = args.inputString
    all_files = glob.glob(os.path.join(resultDir, "*%s" % (inputString)))     # advisable to use os.path.join as this makes concatenation OS independent
    df = pandas.DataFrame()
    
    for f in all_files:
        d = pandas.read_csv(f, comment = "#", sep = ":", header = runs+1)
        d.loc[:,'feature'] = os.path.basename(str(f))
        df  = df.append(d, ignore_index=True)

    df.columns = ["columns","values","feature"]

    df.loc[:,"columns"] = df.loc[:,"columns"].str.replace("\t","").str.replace(" ", "_")

    d = df.pivot(columns = "columns", values = "values", index = "feature").reset_index()

    
    d.loc[:,'feature'].replace(inputString, "", inplace = True, regex = True)

    d.loc[:,'cell'], d.loc[:,'annotation'], d.loc[:,'feature'], d.loc[:,'runParam'] = d.loc[:,'feature'].str.split('.').str

    d.to_csv(args.outputfile, sep="\t", index=False, na_rep="NA")



