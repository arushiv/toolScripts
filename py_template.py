# Standard python template for opening files 
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

if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('mapfile', help="""The file with hg19 tiles to search for.""")
    parser.add_argument('mapchrom', help="""The chromosome for data in the mapfile, this script uses mapfile for one chromosome for speed purposes""")
    parser.add_argument('segmentfile', help="""The file with presence of strong enhancers at 200 bp tile is to be checked.""")
    parser.add_argument('outputfile', help="""Name of the output file corresponding to each mapfile""")
    args = parser.parse_args()
    chrom = str(args.mapchrom)


    with open_maybe_gzipped(args.mapfile) as m_file, open_maybe_gzipped(args.segmentfile) as s_file, open(args.outputfile, 'w+') as o_file:
        map_file = csv.reader(m_file, delimiter='\t', dialect='excel-tab')
        main_list = []
        next(segment_file, None) # skip track header
        for line in segment_file:
           segment_chrom[line[0]].append([line[0], float(line[1]), float(line[2])) # Divide segmentation data by chromosome
        for filename in glob.glob('strongEnhancerSegments_allCells/*.bed'):
                with open_maybe_gzipped(filename) as s_file:
                segment_file = csv.reader(s_file, delimiter='\t', dialect='excel-tab')
                segment_chrom = collections.defaultdict(list)
                if line_start >= segment_start and line_end <= segment_end:    # Falls fully inside
                    fraction += (line_end - line_start) / 200
                    print "found"
                    print fraction
                    break
                elif (segment_start <= line_start <= segment_end) and (line_end > segment_end):   # end coordinate outside segment end coordinate
                    fraction += (segment_end - line_start) / 200
                    line_start = segment_end + 1
                    print "here"
                    print fraction
                elif (line_start < segment_start) and (segment_start <= line_end <= segment_end):   # start coordinate outside segment start coordinate
                    fraction += (line_end - segment_start) / 200
                    print "or here"
                    print fraction
                    break
        if fraction < 0.5:
            temp_list.append(0)
        elif fraction >= 0.5:
            temp_list.append(1)
                
