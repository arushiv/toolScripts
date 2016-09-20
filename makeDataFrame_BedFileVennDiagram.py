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
        return (d['end'] - d['start']).sum(axis=0)
        # return "%.2f" % round(value,2)
    
def intersect2(a, b):
    return a.intersect(b)

def intersect3(a, b, c):
    return a.intersect(b.intersect(c))

def intersect4(a, b, c, d):
    return a.intersect(b.intersect(c.intersect(d)))

def intersect5(a, b, c, d, e):
    return a.intersect(b.intersect(c.intersect(d.intersect(e))))

def form(x):
    string = re.sub('[()]', '', str(x))
    # string = re.sub(r'\)', '', string)
    string = re.sub(',', '', string)
    string = re.sub(' ', '', string)
    return "%s" % (string)

def formBed(x):
    string = re.sub('[()<>]', '', str(x))
    string = re.sub(' ', '', string)
    string = re.sub('BedTool', '', string)
    return "%s" % (string)
    
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="""Script to calculate length of intersecting region for all combinations of up to 5 bed files.""", usage="~/myEnv/bin/python makeDataFrame.py GM12878.broadDomains.annotations.bed,GM12878.hotRegions.annotations.bed,GM12878.stretchEnhancer.annotations.bed,GM12878.superEnhancers.annotations.bed,GM12878.typicalEnhancers.annotations.bed")
    parser.add_argument('bedfile', type=str, help="""bedfile paths, comma separated""")
    args = parser.parse_args()

    
    bedfile = (args.bedfile).split(",")
    
    namelist = [ fetchNames(x) for x in bedfile ]

    beds = [ makeBedTool(x) for x in bedfile ]  # Bedtools for each file

    stuff = list(range(1, len(bedfile)))

    venn_arguments = []
    genome = 3095677412
    
    for a,number in zip(range(0, len(beds)+1), range(0, len(stuff)+1)):
        for subset, subset_number in zip(itertools.combinations(beds, a), itertools.combinations(stuff, number)):
            if len(subset) == 1:
                # print "area%s=%s\t%s" % (form(subset_number), l(*subset), formBed(subset))
                venn_arguments.append("area%s=%s" % (form(subset_number), l(*subset)))

            elif len(subset) == 2:
                # print "n%s=%s\t%s" % (form(subset_number), l(intersect2(*subset)), formBed(subset))
                venn_arguments.append("n%s=%s"  % (form(subset_number), l(intersect2(*subset))))
                
            elif len(subset) == 3:
                # print "n%s=%s\t%s" % (form(subset_number), l(intersect3(*subset)), formBed(subset))
                venn_arguments.append("n%s=%s" % (form(subset_number), l(intersect3(*subset))))
                
            elif len(subset) == 4:
                # print "n%s=%s\t%s" % (form(subset_number), l(intersect4(*subset)), formBed(subset))
                venn_arguments.append("n%s=%s" % (form(subset_number), l(intersect4(*subset))))
                
            elif len(subset) == 5:
                # print "n%s=%s\t%s" % (form(subset_number), l(intersect5(*subset)), formBed(subset))
                venn_arguments.append("n%s=%s" % (form(subset_number), l(intersect5(*subset))))

    print '%s' % ', '.join(map(str, venn_arguments))
    print '%s' % ', '.join(map(str, namelist))
    

    # print bedfile
    # print type(bedfile)
    # pybedtools.contrib.venn_maker.venn_maker(bedfile, names=namelist, figure_filename="test.tiff", script_filename="plot.vennDiagram.R", run=False)
        
