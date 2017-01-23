# Script function:
# 1. For a given matrix of say, RPKM values for different GTEx samples, calculate mean for samples of each tissue. Samples for each tissue specified in another file
# 2. Using the mean RPKM matrix, compute ESI for each tissue

#!/usr/bin/env python
import argparse
import pandas
import math
import glob
import os
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

def log_with_nan(x, y):
    try:
        return math.log(x, y)
    except ValueError:         ## Places with log(0)
        return float('nan')

path = r'.'                     # use your path
all_files = glob.glob("*.dat.gz")     # advisable to use os.path.join as this makes concatenation OS independent

df = pandas.DataFrame()

for f in all_files:
    d = pandas.read_csv(f, comment="#", sep="\t")
    d.loc[:,'feature'] = str(f)
    df  = df.append(d, ignore_index=True)

    
df = df.loc[:,['feature','observed','expected', 'stddev', 'qvalue', 'pvalue']]    
df.loc[:,'feature'].replace(".gatres.dat.gz", "", inplace = True, regex = True)
df.loc[:,'feature'].replace(".annotations", "", inplace = True, regex = True)

df.loc[:,'cell'], df.loc[:,'annotation'], df.loc[:,'feature'] = df.loc[:,'feature'].str.split('.').str

df.loc[:,'enrichment'] = (df.loc[:,'observed'] / df.loc[:,'expected']).apply(lambda x: log_with_nan(x, 2))

df.to_csv("results.dat", sep="\t", index=False, na_rep="NA")
p = sns.FacetGrid(df, col="cell",  hue="annotation").map(sns.stripplot, "feature", "enrichment")

p.savefig("fig.pdf")
