# Comment 1: For each 200bp tile of hg19 we have lineinates in mouse or rat by running bnMapper, filtered for tiles which map both to mouse and rat. Use this script to fetch the chromatin state segmentatations corresponding to these lineinates in particular mouse and rat tissue sample datasets
# Comment 2: Have confirmed earlier using checkMappingRegion.py script that all hg19 tiles map to more or less contiguous regions in mouse and rat; ie chromosome is same for all sub mappings of a particular tile
# The map file is of the form: human_tiles; #_mappings_mouse; total_length_mappings_mouse; #_mappings_rat; total_length_mappings_rat; lineinates_mouse_comma_separated; lineinates_rat_comma_separated
# The segment file should be ascending order sorted
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

    parser.add_argument('mapfile', help="""The file with hg19 tiles and mouse and rat lineinates they map to.""")
    parser.add_argument('segmentfile', help="""The file with chromatin segmentations to be fetched.""")
    parser.add_argument('organism', help="""Organism for cell/tissue segmentation file""")
    args = parser.parse_args()

    org = str(args.organism)
    with open_maybe_gzipped(args.mapfile) as m_file, open_maybe_gzipped(args.segmentfile) as s_file:

        map_file = csv.reader(m_file, delimiter=';', dialect='excel-tab')
        segment_file = csv.reader(s_file, delimiter='\t', dialect='excel-tab')
        segment_chrom = collections.defaultdict(list)


        for line in segment_file:
            segment_chrom[line[0]].append([line[0], float(line[1]), float(line[2]), line[3]])

        if str(org) == 'mouse':
            for line in map_file:
                fraction = collections.defaultdict(float)
                start = 6
                end = 7
                minisegments_length = float(line[2])
                for nc in range(int(line[1])):             # While we check every submapping for each tile
                    line_start = float(line[start])
                    line_end = float(line[end])
                    for segment in segment_chrom[line[5]]:    # Refer to comment 2
                        segment_start = segment[1]
                        segment_end = segment[2]
                        segment_key = segment[3]
                        if line_start >= segment_start and line_end <= segment_end:    # Falls fully inside
                            fraction[segment_key] += (line_end - line_start) / minisegments_length
                            start += 3
                            end += 3
                            break
                        elif (segment_start <= line_start <= segment_end) and (line_end > segment_end):   # end coordinate outside segment end coordinate
                            fraction[segment_key] += (segment_end - line_start) / minisegments_length
                            line_start = segment_end + 1
                        elif (line_start < segment_start) and (segment_start <= line_end <= segment_end):   # start coordinate outside segment start coordinate
                            fraction[segment_key] += (line_end - segment_start) / minisegments_length
                            start += 3
                            end += 3
                            break
                if not fraction:
                    print(' {} {} {} not_found'.format(line[0], line[1], line[2]))
                else:
                    print(line[0], line[1], line[2], fraction, max(fraction.iteritems(), key=operator.itemgetter(1)))

        if str(org) == 'rat':
            for line in map_file:
                fraction = collections.defaultdict(float)
                chrom = 5 + (int(line[1])*3)
                start = 6 + (int(line[1])*3)
                end = start + 1
                minisegments_length = float(line[4])
                for nc in range(int(line[3])):             # While we check every submapping for each tile
                    line_start = float(line[start])
                    line_end = float(line[end])
                    for segment in segment_chrom[line[chrom]]:    # Refer to comment 2
                        segment_start = segment[1]
                        segment_end = segment[2]
                        segment_key = segment[3]
                        if line_start >= segment_start and line_end <= segment_end:    # Falls fully inside
                            fraction[segment_key] += (line_end - line_start) / minisegments_length
                            start += 3
                            end += 3
                            break
                        elif (segment_start <= line_start <= segment_end) and (line_end > segment_end):   # end coordinate outside segment end coordinate
                            fraction[segment_key] += (segment_end - line_start) / minisegments_length
                            line_start = segment_end + 1
                        elif (line_start < segment_start) and (segment_start <= line_end <= segment_end):   # start coordinate outside segment start coordinate
                            fraction[segment_key] += (line_end - segment_start) / minisegments_length
                            start += 3
                            end += 3
                            break
                if not fraction:
                    print('{} {} {} not_found'.format(line[0], line[3], line[4]))
                else:
                    print(line[0], line[3], line[4], fraction, max(fraction.iteritems(), key=operator.itemgetter(1)))
