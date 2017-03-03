# Script function:
# 1. For a given matrix of say, RPKM values for different GTEx samples, calculate mean for samples of each tissue. Samples for each tissue specified in another file


#!/usr/bin/env python
import argparse
import pandas
import math


def log_with_nan(x, y):
    try:
        return math.log(x, y)
    except ValueError:         ## Places with log(0)
        return float('nan')

def entropy_function(x):
    y = -x * log_with_nan(x, 2)
    return y


def calculate_ESI(df):
    d = df.loc[:,["mean.Pituitary","mean.Spleen","mean.Brain.CerebellarHemisphere","mean.Prostate","mean.Brain.FrontalCortex.BA9.","mean.Brain.Nucleusaccumbens.basalganglia.","mean.Brain.Cortex","mean.Brain.Caudate.basalganglia.","mean.Cells.EBV.transformedlymphocytes","mean.Liver","mean.Brain.Cerebellum","mean.Artery.Coronary","mean.AdrenalGland","mean.Colon.Sigmoid","mean.Esophagus.GastroesophagealJunction","mean.Pancreas","mean.Testis","mean.Stomach","mean.Heart.AtrialAppendage","mean.Colon.Transverse","mean.Breast.MammaryTissue","mean.Heart.LeftVentricle","mean.Artery.Aorta","mean.Adipose.Visceral.Omentum.","mean.Esophagus.Muscularis","mean.Skin.NotSunExposed.Suprapubic.","mean.Cells.Transformedfibroblasts","mean.Esophagus.Mucosa","mean.Nerve.Tibial","mean.Lung","mean.Thyroid","mean.Artery.Tibial","mean.Adipose.Subcutaneous","mean.Skin.SunExposed.Lowerleg.","mean.WholeBlood","mean.Muscle.Skeletal"]]
    d_relative_rpkm = d.div(d.sum(axis=1), axis=0)  ## Relative RPKM = x(g,t)/sum(x(g,t))
    d.loc[:,'entropy'] = d_relative_rpkm.applymap(lambda x: entropy_function(x)).sum(axis=1)    ## Entropy = - sum(relativeRPKM(g,t) * log2(relative RPKM(g,t)))
    d_qvalues = d_relative_rpkm.applymap(lambda x: log_with_nan(x, 2)).apply(lambda x: d.loc[:,'entropy'] - x, axis=0)   ## Q(g,t) = Entropy(g) - log2(relativeRPKM(g,t))
    d_max_q = d_qvalues.max(axis=0).to_frame().transpose()        ## maxQ = max(Q(t))

    d_esi = d_qvalues.apply(lambda x: 1 - x/d_max_q.squeeze(), axis=1)
    return d_esi

    
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='For a given matrix of say, RPKM values for different GTEx samples, calculate mean for samples of each tissue. Samples for each tissue specified in another file. Then calculate ESI for each tissue', usage='python esiScoreAfterMean.py datafile.txt samplefile.txt output.txt')
    parser.add_argument('datafile', help="""Input Matrix. Tab separated, with header specifying sample ID. First two columns taken as 'gene name' and 'description'. Rest columns as different sample numbers""")
    # parser.add_argument('samplefile', help="""Tab separated file with 2 columns: 'tissue' 'sample', specifying sample IDs corresponding to each tissue type.""")
    parser.add_argument('outputfile', help ="""Output file.""")
    args = parser.parse_args()

    datafile = pandas.read_csv(args.datafile, sep='\t', index_col=[0,1])
    # samplefile = pandas.read_csv(args.samplefile, sep='\t')
    outputfile = args.outputfile

    # Calculate Mean for each tissue given the sample names 
    # sampledf = samplefile.groupby('tissue')['sample'].apply(list)
    # df = sampledf.apply(lambda x: datafile.loc[:,x].mean(axis=1)).transpose()#.reset_index()

    # Calculate ESI
    d_esi = calculate_ESI(datafile)
    d_esi.to_csv(outputfile, sep='\t', na_rep="NA")
