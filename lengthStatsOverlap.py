#version of ./unusedSciptsAndOutputs/lengthStatsAndOverlap.py optimized by John Hensley
## For each 200bp tile of hg19 (hg19_200tile.bed4.gz), search mappings in maps to mm9 (hg19tileMapTomm9.bed) and rn5 (hg19tileMapTorn5.bed), calculate number of maps in each, total length of maps and coordinates of individual mappings
#!/usr/bin/env python

from __future__ import print_function
import argparse
import collections
import csv
import gzip
import magic
import pprint

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

    parser.add_argument('human', help="""The human file.""")
    parser.add_argument('mouse', help="""The mouse file.""")
    parser.add_argument('rat', help="""The rat file.""")

    args = parser.parse_args()

    with open_maybe_gzipped(args.human) as h_file, open_maybe_gzipped(args.mouse) as m_file, open_maybe_gzipped(args.rat) as r_file:
        human_file = csv.reader(h_file, dialect='excel-tab')
        mouse_file = csv.reader(m_file, dialect='excel-tab')
        rat_file = csv.reader(r_file, dialect='excel-tab')

        mouse = collections.defaultdict(list)
        rat = collections.defaultdict(list)

        for line in mouse_file:
            mouse[line[3]].append(line)

        for line in rat_file:
            rat[line[3]].append(line)

        for line in human_file:
            key = line[3]
            mouse_matches = mouse[key]
            mouse_total_len = sum([abs(int(mouse_line[2]) - int(mouse_line[1])) for mouse_line in mouse_matches])
            mouse_coords = [','.join(mouse_line[0:3]) for mouse_line in mouse_matches]

            rat_matches = rat[key]
            rat_total_len = sum([abs(int(rat_line[2]) - int(rat_line[1])) for rat_line in rat_matches])
            rat_coords = [','.join(rat_line[0:3]) for rat_line in rat_matches]

            print('; '.join(map(str, [key, len(mouse_matches), mouse_total_len, len(rat_matches), rat_total_len, mouse_coords, rat_coords])))
