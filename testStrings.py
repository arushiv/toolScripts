import sys
import fileinput

for line in fileinput.input():
	print repr(str(line)) 
