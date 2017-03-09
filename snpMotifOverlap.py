# Find overlap position and probabilities and information content at that position  
#!/usr/bin/env python
# snp_chrom snp_start snp_end motif_chrom motif_start motif_end motif strand
# from __future__ import print_function
import argparse
import collections
import csv
import gzip
import magic
import pprint
import re
import operator
import glob
import math
import sys
import subprocess as sp

def open_maybe_gzipped(f):
    with open(f) as test_read:
        mimetype = magic.from_buffer(test_read.read(1024), mime=True).decode('utf-8')
    if mimetype == 'application/x-gzip':
        f = gzip.open(f, mode='rt') # Python3: , encoding='utf-8')
    else:
        f = open(f, 'rt')
    return f

def calculate_overlap_position(x,y):
    return x-y

def calculate_overlap_position_forward_logo(a,b,c,d):
    if d == "+" :
        return a-b
    elif d == "-":
        return c - a + 1

def calculate_motif_probabilites(x,y,z):
    with open("/home/porchard/mats/" + x + ".meme") as motif_file:      # Path for motif pwm .meme files
        for n,line in enumerate(motif_file):
            if "probability" in line:
                zeroth = n
        motif_file.seek(0)
        for i,line in enumerate(motif_file):
            if i == zeroth+y:
                probability = line
        if z=="+":
            return probability.rstrip()
        elif z=="-":
            p = probability.split()
            return p[3]+'\t'+p[2]+'\t'+p[1]+'\t'+p[0]

def calculate_motif_information_content(x):
    x=x.split()
    x = ",".join(x)
    cmd = "echo %s | perl /home/scjp/src/parker/tools/informationColumnCommaDelimited.pl -f - -p" %x     # Path for Steve's script to calculate Information content at given position
    return float((sp.check_output(cmd, shell = True)).split()[-1])*2

def checkFunc(motif, overlap_pos_forwardlogo, strand):
    try:
        return calculate_motif_probabilites(motif, overlap_pos_forwardlogo, strand)
    except UnboundLocalError:         ## Places with log(0)
        return str(0)+'\t'+str(0)+'\t'+str(0)+'\t'+str(0)
    
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Take SNP coordinates and TF footprint coordinates along with footprint strand and motif PWM information, report the probabilites to A T G and C bases in the motif at the SNP-Motif overlap position. Also output information content at the position. Adjust path to motif pwm .meme files in the script.', usage='python snpMotifOverlap.py <datafile>' )
    parser.add_argument('datfile', help="""The data file, tab delimited. Should contain at least these columns:\nsnp_chrom snp_start snp_end motif_chrom motif_start motif_end motif strand. \n Output has the following columns added to the input file: overlap_position overlap_position_relative_to_fo\
rward_logo motif_prob_A motif_prob_C    motif_prob_G    motif_prob_T    pwinformation_content""")
    args = parser.parse_args()


    
    with open_maybe_gzipped(args.datfile) as d_file:
        dat_file = csv.reader(d_file, delimiter='\t', dialect='excel-tab')
        for line in dat_file:
            snp_pos = int(line[2])
            motif_start = int(line[4])
            motif_end = int(line[5])
            motif = str(line[6])
            strand = str(line[7])
            overlap_pos = calculate_overlap_position(snp_pos,motif_start)
            overlap_pos_forwardlogo = calculate_overlap_position_forward_logo(snp_pos, motif_start, motif_end, strand)
            
            motif_probabilities = checkFunc(motif, overlap_pos_forwardlogo, strand)
            # motif_probabilities = calculate_motif_probabilites(motif, overlap_pos_forwardlogo, strand)
            # print motif_probabilities
            motif_information_content = calculate_motif_information_content(motif_probabilities)
            # print motif_information_content.split()[-1]
            line.extend([str(overlap_pos), str(overlap_pos_forwardlogo), str(motif_probabilities), str(motif_information_content)])
            print '\t'.join(line)

