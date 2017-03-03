# Calculates how many mappings map to mouse and rat.
import sys
import re
import gzip
import fileinput
from collections import defaultdict
from pylab import *
import matplotlib.pyplot as plt

ch_map_total = defaultdict(int)
ch_map_mouseonly = defaultdict(int)
ch_map_ratonly = defaultdict(int)
ch_map_both = defaultdict(int)
ch_map_none = defaultdict(int)
#map_len_mouse = []
#map_len_rat = []

print 'human_chromosome\ttotal_tiles\tmouse_only_mappings\tpercent_mouse_only_mappings\trat_only_mappings\tpercent_rat_only_mappings\tmappings_mouseAndrat\tpercent_mappings_mouseAndrat\tunmapped_tiles\tpercent_unmapped_tiles'
for l in fileinput.input():
    lt = l.rstrip()
    x = re.split(r'[":"";"]',lt)
    ch_map_total[str(x[0])] += 1
 #   map_len_mouse.append(int(x[4]))
  #  map_len_rat.append(int(x[6]))
    
    if int(x[3]) > 0 and int(x[5]) == 0:
        #print 'Exact %s' % (l)
        ch_map_mouseonly[str(x[0])] += 1
#        print ch_map_mouseonly[str(x[0])]
    elif int(x[3]) == 0 and int(x[5]) > 0:
        #print 'Inside %s' % (l)
        ch_map_ratonly[str(x[0])] += 1
 #       print ch_map_ratonly[str(x[0])]
    elif int(x[3]) > 0 and int(x[5]) > 0:
        #print 'Inside/match %s' % (l)
        ch_map_both[str(x[0])] += 1
    elif int(x[3]) == 0 and int(x[5]) == 0:
        ch_map_none[str(x[0])] += 1

            
for m, n in ch_map_total.iteritems():
    print '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (m, n, ch_map_mouseonly[m],(float(ch_map_mouseonly[m])*100/n), ch_map_ratonly[m], (float(ch_map_ratonly[m])*100/n), ch_map_both[m], (float(ch_map_both[m])*100/n), ch_map_none[m], (float(ch_map_none[m])*100/n))
  

#plt.hist(map_len_mouse)
#plt.save('hist_mouse.png')

#plt.hist(map_len_rat)
#plt.save('hist_rat.png')
    
