# Run GREGOR: Fix conf files, submit drmr commands, compile dataframe and plot 
#!/usr/bin/env python
# 
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
import os
import itertools
import shutil
import fileinput
import fnmatch

def stringReplace(filename, findstring, replacestring):
    with open("temp.txt", 'w') as fout:
        with open(filename, 'r') as fin:
            for line in fin:
                if findstring in line:
                    fout.write('%s ## %s'%(replacestring, line.split("##")[-1]))
                else:
                    fout.write(line)
    os.rename("temp.txt", filename)

class RunGregor(object):
    def __init__(self, workingDirectory):
        self.workingDirectory = workingDirectory
    
    def makeOutputDirectories(self, snpFileDirectory, bedFileDirectory, bedFilesPerJob, linkageDisequilibrium, nameGregorRun):
        snpFileCounter = len(glob.glob1(snpFileDirectory,"*.txt"))
        ldCounter=len(linkageDisequilibrium)
        self.nameGregorRun = nameGregorRun
        if bedFilesPerJob==1:
            self.bedCounter=1
        else:
            self.bedCounter=((len(glob.glob1(bedFileDirectory,"*.bed")))/bedFilesPerJob) + 1
        numberOfGregorJobs=snpFileCounter*ldCounter*self.bedCounter

        self.fileList = []
        for snpFile in glob.glob1(snpFileDirectory,"*.txt"):
            for ldValue in linkageDisequilibrium:
                name = snpFile.replace(".txt","")
                if len(linkageDisequilibrium) != 1:
                    name1 = "%s_%s"%(name,ldValue)
                else:
                    name1 = name
                for bedValue in range(self.bedCounter):
                    if self.bedCounter != 1:
                        filename = "%s_%s"%(name1,bedValue+1)
                    else:
                        filename = name1
                    if self.nameGregorRun != "":
                        filename = "%s_%s"%(filename, self.nameGregorRun)
                    newpath = os.path.join(self.workingDirectory, "output_%s"%filename)
                    self.fileList.append(filename)
                    if not os.path.exists(newpath):
                        os.makedirs(newpath)
                    else:
                        print "Caution: Directory already exists!!"
    
    def makeIndexBedFiles(self, bedFileDirectory, bedFilesPerJob):
        # elif bedFilesPerJob=="name":
        #     # find number of unique names
        tempfile=self.workingDirectory + 'temp.txt'
        with open(tempfile,'w') as f:
            for bedFile in glob.glob1(bedFileDirectory,"*.bed"):
                f.write(os.path.join(bedFileDirectory, bedFile)+"\n")
        if self.bedCounter == 1:
            os.rename(tempfile, os.path.join(self.workingDirectory, 'bedFileIndex_%s.txt'%(self.nameGregorRun)))
        else:
            input = open(tempfile,'r').read().split('\n')
            at = 1
            for lines in range(0, len(input), bedFilesPerJob):
                outputData = input[lines:lines+bedFilesPerJob]
                output = open(os.path.join(self.workingDirectory,'bedFileIndex_%s_%s.txt'%(str(at),self.nameGregorRun)), 'w')
                output.write('\n'.join(outputData))
                output.close()
                at += 1
           
    def makeConfFiles(self, snpFileDirectory, linkageDisequilibrium, sampleConfFile, cores):
        for filename in self.fileList:
            namelist = filename.split('_')
            newfile = os.path.join(self.workingDirectory, "enrich_%s.conf"%filename)
            shutil.copy(sampleConfFile, newfile)
            stringReplace(newfile, "INDEX_SNP_FILE =","INDEX_SNP_FILE = %s"%os.path.join(snpFileDirectory,"%s.txt"%namelist[0]))
            if self.bedCounter == 1:
                if self.nameGregorRun != "":
                    stringReplace(newfile, "BED_FILE_INDEX =","BED_FILE_INDEX = %s"%os.path.join(self.workingDirectory,"bedFileIndex_%s.txt"%self.nameGregorRun))
                else:
                    stringReplace(newfile, "BED_FILE_INDEX =","BED_FILE_INDEX = %s"%os.path.join(self.workingDirectory,"bedFileIndex.txt"))
            else:
                if self.nameGregorRun == "":
                    stringReplace(newfile, "BED_FILE_INDEX =","BED_FILE_INDEX = %s"%os.path.join(self.workingDirectory,"bedFileIndex_%s.txt"%namelist[-1]))
                else:
                    stringReplace(newfile, "BED_FILE_INDEX =","BED_FILE_INDEX = %s"%os.path.join(self.workingDirectory,"bedFileIndex_%s_%s.txt"%(namelist[-2],self.nameGregorRun)))
                    
            if len(linkageDisequilibrium) != 1:
                stringReplace(newfile, "R2THRESHOLD =","R2THRESHOLD = %s"%namelist[1])
            stringReplace(newfile, "OUT_DIR =","OUT_DIR = %s"%os.path.join(self.workingDirectory,"output_%s"%filename))
            stringReplace(newfile, "JOBNUMBER =","JOBNUMBER = %s"%cores)            
                         
                  
    def makeDrmrFile(self, cores, walltime):
        self.runfile = os.path.join(self.workingDirectory, 'jobs_%s.dra'%(self.nameGregorRun)) 
        with open(self.runfile, 'w') as f:
            f.write("# drmr:job nodes=1 processors=%s processor_memory=4000 time_limit=%s\n"%(cores, walltime))
            for filename in self.fileList:
                f.write("GREGOR.pl --conf enrich_%s.conf\n"%filename)
            f.write("\n# drmr:wait\n")
            f.write("\n# drmr:job nodes=1 processors=1 processor_memory=4000 time_limit=1:00:00\n")
            for filename in self.fileList:
                f.write(" python /home/arushiv/toolScripts/makeDataFrame_gregor.py output_%s/StatisticSummaryFile.txt stats_%s.dat; "%(filename, self.nameGregorRun))

    def submitDrmrJob(self):
        # runfile = os.path.join(self.workingDirectory, runfile)
        cmd = "drmr %s -j job"%self.runfile
        sp.call(cmd, shell=True)
    
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Run GREGOR: Fix conf files, submit drmr commands, compile dataframe and plot.')
    parser.add_argument('snpFileDirectory', type=str, help="""The directory with SNP files where GREGOR would be run for each such file. Should be in chr:pos or rsID format; filename should not contain '_'""")
    parser.add_argument('bedFileDirectory', type=str, help="""The directory with (unzipped) bed files""")
    parser.add_argument('-f','--sampleConfFile', nargs='?', default='/home/arushiv/toolScripts/gregor_sampleConfFile.conf',
                        help="""Sample .conf file where only snp file, bed file and output folder will be changed. (default: /lab/arushiv/toolScripts/gregor_sampleConfFile.conf)""")
    parser.add_argument('-n','--bedFilesPerJob', nargs='?', default='1', help="""Number of bedfiles per GREGOR job. (default: GREGOR runs all bedfiles in one job for each SNP file)""")
    parser.add_argument('-ld', '--linkageDisequilibrium', nargs='?', default='0.8', help="""minimum LD for proxy SNPs. Multiple LD values to be supplied comma (,) separated. (default: 0.8)""")
    parser.add_argument('-t','--walltime', nargs='?', type=str, default='12:00:00', help="""Walltime for each GREGOR job. (default: 12:00:00 [hh:mm:ss])""")
    parser.add_argument('-s','--cores', nargs='?', type=str, default='6', help="""Number of cores for each GREGOR job. (default: 6)""")
    parser.add_argument('-d','--workingDirectory', type=str, nargs='?', default='.', help="""Directory where output directories will be stored. (default: current directory)""")
    parser.add_argument('--now', action='store_const', const='now', default='wait', help="""Submit drmr job file now.""")
    parser.add_argument('-name','--nameGregorRun', type=str, nargs='?', default="", help="""Name string to be appended to conf files and output folders indicating specific gregor run. (default: "")""")
    args = parser.parse_args()

    snpFileDirectory=args.snpFileDirectory
    bedFileDirectory=args.bedFileDirectory
    sampleConfFile=args.sampleConfFile
    if args.bedFilesPerJob != "name":
        bedFilesPerJob=int(args.bedFilesPerJob)
    else:
        bedFilesPerJob="name"
    ld = (args.linkageDisequilibrium).split(',')
    walltime = args.walltime
    cores = args.cores
    workingDirectory = args.workingDirectory
    runnow = args.now
    nameGregorRun = args.nameGregorRun

    f = RunGregor(workingDirectory)
    f.makeOutputDirectories(snpFileDirectory, bedFileDirectory, bedFilesPerJob, ld, nameGregorRun)
    f.makeIndexBedFiles(bedFileDirectory, bedFilesPerJob)
    f.makeConfFiles(snpFileDirectory, ld, sampleConfFile, cores)
    f.makeDrmrFile(cores, walltime)

    if runnow == "now":
        f.submitDrmrJob()
    
    # makeDrmrFile()


