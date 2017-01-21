# Script function: For given bed files (up to 5), print all combinations of intersection lengths
# Output can be used to construct venn diagrams

#!/usr/bin/env python
import argparse
import pybedtools
import glob
import itertools
import re
import os
import pandas
import subprocess as sp

def fetchNames(x):
    string = os.path.splitext(os.path.basename(x))[0]
    return string.split('.')[1]

def makeBedTool(x):
    return pybedtools.BedTool(x)                # Columns in input dataframe 'chr,  start,  end, name, score, strand,  motif'


def l(x):
    if x.count() == 0:
        return 0
    else:
        d = x.to_dataframe()
        return float((d['end'] - d['start']).sum(axis=0))
        # return "%.2f" % round(value,2)
    
def intersect2(a, b):
    return a.intersect(b)

def intersect3(a, b, c):
    return a.intersect(b.intersect(c))

def intersect4(a, b, c, d):
    return a.intersect(b.intersect(c.intersect(d)))

def intersect5(a, b, c, d, e):
    return a.intersect(b.intersect(c.intersect(d.intersect(e))))

def form(x,n):
    string = re.sub('[()]', '', str(x))
    string = re.sub("'", "", string)
    string = re.sub(',', '', string)
    string = re.sub(' ', '&', string)
    if n == 1:
        return "%s" % (string)
    else:
        return "`%s`" % (string)

def makeRFile(printString, Rfile):
    Rfile=open(Rfile, 'w+')
    Rfile.write("library(UpSetR)\n\n")
    Rfile.write("expressionInput <- c(%s)\n\n" %(printString))
    Rfile.write("pdf('fig.pdf', width=4, height=4)\n")
    Rfile.write("upset(fromExpression(expressionInput), order.by='freq')\n")
    Rfile.write("dev.off()")
        
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="""Script to calculate length of intersecting region for all combinations of up to 5 bed files.""", usage="~/myEnv/bin/python makeDataFrame.py GM12878.broadDomains.annotations.bed,GM12878.hotRegions.annotations.bed,GM12878.stretchEnhancer.annotations.bed,GM12878.superEnhancers.annotations.bed,GM12878.typicalEnhancers.annotations.bed")
    parser.add_argument('bedfile', type=str, help="""bedfile paths, comma separated""")
    parser.add_argument('Rfile', type=str, help="""R script to run.""")
    parser.add_argument('--now', action='store_const', const='now', default='wait', help="""Run R script now. It saves a 4*4 figure at fig.pdf, Plot ordered by frequency.""")
    args = parser.parse_args()
    
    bedfile = (args.bedfile).split(",")
    Rfile = args.Rfile
    
    namelist = [ fetchNames(x) for x in bedfile ]
    beds = [ makeBedTool(x) for x in bedfile ]  # Bedtools for list(range(1, len(bedfile)+1))

    stuff = list(range(1, len(bedfile)+1))

    venn_arguments = []

    for a,number in zip(range(0, len(beds)+1), range(0, len(stuff)+1)):
        for subset, subset_name in zip(itertools.combinations(beds, a), itertools.combinations(namelist, number)):
          
            if len(subset) == 1:
                venn_arguments.append("%s=%s" % (form(subset_name, 1), l(*subset)/1000000))
                
            elif len(subset) == 2:
                venn_arguments.append("%s=%s"  % (form(subset_name, 2), l(intersect2(*subset))/1000000))
                
            elif len(subset) == 3:
                venn_arguments.append("%s=%s" % (form(subset_name, 3), l(intersect3(*subset))/1000000))
                
            elif len(subset) == 4:
                venn_arguments.append("%s=%s" % (form(subset_name, 4), l(intersect4(*subset))/1000000))
                
            elif len(subset) == 5:
                venn_arguments.append("%s=%s" % (form(subset_name, 5), l(intersect5(*subset))/1000000))

    printString = '%s' % ', '.join(map(str, venn_arguments))

    makeRFile(printString, Rfile)
    print "R script created at %s" %(Rfile)
    

    if args.now == "now":
        cmd = "Rscript %s" %(Rfile)
        sp.call(cmd, shell=True)
        print "Figure saved at fig.pdf"
