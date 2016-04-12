# Match gene name on provided column, separate into files as per quantiles.
# Ex. binQuintile file at /lab/arushiv/chromatin/2015_12_08_islet_eQTL/eqtl_binnedByIsletGenes/bins.protein_coding.binQuintile.1.byGenes contains gene names that are in Quintile 1
#!/usr/bin/env python


#from __future__ import print_function
import argparse
import collections
import csv
import gzip
import magic
import pprint
import re
import operator
import glob

def open_maybe_gzipped(f):
    with open(f) as test_read:
        mimetype = magic.from_buffer(test_read.read(1024), mime=True).decode('utf-8')
    if mimetype == 'application/x-gzip':
        f = gzip.open(f, mode='rt') # Python3: , encoding='utf-8')
    else:
        f = open(f, 'rt')
    return f

def splitByQuantile(quantilef, eqtlf, column_gene, qf, outputfile):
    with open_maybe_gzipped(eqtlf) as efile:
        eqtlfile = csv.reader(efile, delimiter='\t', dialect='excel-tab')

        with open_maybe_gzipped(quantilef) as qfile:
            quantilefile = csv.reader(qfile, delimiter='\t', dialect='excel-tab')

            with open(outputfile, 'w') as fout:
                for line in eqtlfile:
                    for stuff in quantilefile:
                        if (line[column_gene - 1]).split('.')[0] == stuff[qf - 1]:
                            fout.write('\t'.join(line))
                            fout.write('\n')
                    qfile.seek(0)
                    
if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('quantilefile', help="""The file with genes for a certain quantile.""")
    parser.add_argument('eqtlfile', help="""The file with data to be searched and split.""")
    parser.add_argument('-f', '--column_gene', type=int, default=1, help="""The column number in the eqtl file which contains the gene name.""")
    parser.add_argument('-qf', '--column_quantileFile', type=int, default=1, help="""The column number in quantile file for the gene name.""")
    parser.add_argument('outputfile', help="""Name of the output file.""")

    args = parser.parse_args()

    quantilef = args.quantilefile
    eqtlf = args.eqtlfile
    column_gene = args.column_gene
    outputfile = args.outputfile
    qf = args.column_quantileFile
    
    splitByQuantile(quantilef, eqtlf, column_gene, qf, outputfile)
    
    # with open_maybe_gzipped(args.mapfile) as m_file, open_maybe_gzipped(args.segmentfile) as s_file, open(args.outputfile, 'w+') as o_file:
    #     map_file = csv.reader(m_file, delimiter='\t', dialect='excel-tab')
    #     main_list = []
    #     next(segment_file, None) # skip track header
    #     for line in segment_file:
    #        segment_chrom[line[0]].append([line[0], float(line[1]), float(line[2])) # Divide segmentation data by chromosome
    #     for filename in glob.glob('strongEnhancerSegments_allCells/*.bed'):
    #             with open_maybe_gzipped(filename) as s_file:
    #             segment_file = csv.reader(s_file, delimiter='\t', dialect='excel-tab')
    #             segment_chrom = collections.defaultdict(list)
    #             if line_start >= segment_start and line_end <= segment_end:    # Falls fully inside
    #                 fraction += (line_end - line_start) / 200
    #                 print "found"
    #                 print fraction
    #                 break
    #             elif (segment_start <= line_start <= segment_end) and (line_end > segment_end):   # end coordinate outside segment end coordinate
    #                 fraction += (segment_end - line_start) / 200
    #                 line_start = segment_end + 1
    #                 print "here"
    #                 print fraction
    #             elif (line_start < segment_start) and (segment_start <= line_end <= segment_end):   # start coordinate outside segment start coordinate
    #                 fraction += (line_end - segment_start) / 200
    #                 print "or here"
    #                 print fraction
    #                 break
    #     if fraction < 0.5:
    #         temp_list.append(0)
    #     elif fraction >= 0.5:
    #         temp_list.append(1)
                
