# For each 200bp tile of hg19 we have coordinates in mouse or rat by running bnMapper. Use this script to filter out tiles which map both to mouse and rat. Use another script to fetch the chromatin state segmentatations corresponding to these coordinates in particular mouse and rat tissue sample datasets  
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

    parser.add_argument('humantile', help="""The humantile file.""")
    parser.add_argument('mouse', help="""The mouse file.""")
    parser.add_argument('rat', help="""The rat file.""")
    parser.add_argument('output', help="""The dual mappable output filename""")
    parser.add_argument('totalout', help="""The total output filename""" )
    args = parser.parse_args()

    with open_maybe_gzipped(args.humantile) as h_file, open_maybe_gzipped(args.mouse) as m_file, open_maybe_gzipped(args.rat) as r_file, open(args.output, 'w+') as out_file, open(args.totalout, 'w+') as totalout_file:
        humantile_file = csv.reader(h_file, dialect='excel-tab')
        mouse_file = csv.reader(m_file, dialect='excel-tab')
        rat_file = csv.reader(r_file, dialect='excel-tab')

        mouse = collections.defaultdict(list)
        rat = collections.defaultdict(list)

        for line in mouse_file:
            mouse[line[3]].append(line)

        for line in rat_file:
            rat[line[3]].append(line)

        for line in humantile_file:
            key = line[3]
            mouse_matches = mouse[key]
            mouse_total_len = sum([abs(int(mouse_line[2]) - int(mouse_line[1])) for mouse_line in mouse_matches])
            mouse_coords = [','.join(mouse_line[0:3]) for mouse_line in mouse_matches]

            rat_matches = rat[key]
            rat_total_len = sum([abs(int(rat_line[2]) - int(rat_line[1])) for rat_line in rat_matches])
            rat_coords = [','.join(rat_line[0:3]) for rat_line in rat_matches]

            totalout_file.write('; '.join(map(str, [key, len(mouse_matches), mouse_total_len, len(rat_matches), rat_total_len, mouse_coords, rat_coords, '\n'])))

            if len(mouse_matches) > 0 and len(rat_matches) > 0: 
                out_file.write('; '.join(map(str, [key, len(mouse_matches), mouse_total_len, len(rat_matches), rat_total_len, mouse_coords, rat_coords, '\n'])))
        
            
#print '%s; %s; %s; %s; %s; %s; %s' % (sh[3], n_m, total_len_m, n_r, total_len_r, coords_m, coords_r)
