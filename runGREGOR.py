# Run GREGOR: Fix conf files, submit drmr commands, compile dataframe and plot 
#!/usr/bin/env python


import argparse
import glob
import sys
import subprocess as sp
import os
import shutil


# To find and replace parameter values to create .conf files accroding to supplied parameters
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

# To make GREGOR output directories for each .conf file
    def makeOutputDirectories(self, snpFileDirectory, bedFileDirectory, bedFilesPerJob, linkageDisequilibrium, nameGregorRun):
        snpFileCounter = len(glob.glob1(snpFileDirectory,"*.txt"))
        ldCounter=len(linkageDisequilibrium)
        self.nameGregorRun = nameGregorRun
        if bedFilesPerJob=="all":
            self.bedCounter=1
        elif len(glob.glob1(bedFileDirectory,"*.bed")) % bedFilesPerJob == 0:
            self.bedCounter=((len(glob.glob1(bedFileDirectory,"*.bed")))/bedFilesPerJob)
            print self.bedCounter
        else:
            self.bedCounter=((len(glob.glob1(bedFileDirectory,"*.bed")))/bedFilesPerJob) + 1
            print self.bedCounter
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
            elif len(linkageDisequilibrium) == 1 and linkageDisequilibrium[0] != "0.8":
                stringReplace(newfile, "R2THRESHOLD =","R2THRESHOLD = %s"%linkageDisequilibrium[0])
                        
            stringReplace(newfile, "OUT_DIR =","OUT_DIR = %s"%os.path.join(self.workingDirectory,"output_%s"%filename))
            stringReplace(newfile, "JOBNUMBER =","JOBNUMBER = %s"%cores)            
                         
# Make file to submit jobs in the scheduler                  
    def makeDrmrFile(self, cores, walltime):
        self.runfile = os.path.join(self.workingDirectory, 'jobs_%s.dra'%(self.nameGregorRun)) 
        with open(self.runfile, 'w') as f:
            f.write("# drmr:job nodes=1 processors=%s processor_memory=4000 time_limit=%s\n"%(cores, walltime))  # GREGOR job resource parameters
            for filename in self.fileList:
                f.write("GREGOR.pl --conf enrich_%s.conf\n"%filename)                                            # Print one job command for each .conf file 
            f.write("\n# drmr:wait\n")                                                                           # Wait till all GREGOR jobs finish running
            f.write("\n# drmr:job nodes=1 processors=1 processor_memory=4000 time_limit=1:00:00\n")              # Resource parameters to assemble data frame
            for filename in self.fileList:
                f.write(" python /home/arushiv/toolScripts/makeDataFrame_gregor.py output_%s/StatisticSummaryFile.txt stats_%s.dat; "%(filename, self.nameGregorRun))   # Assemble dataframe from result StatisticSummaryFile.txt in each output folder. Supply path to script makeDataFrame_gregor.py here

    def submitDrmrJob(self):
        cmd = "drmr %s -j job"%self.runfile     # Syntax to submit drmr job
        sp.call(cmd, shell=True)
    
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Run GREGOR: Fix conf files, submit GREGOR jobs using drmr commands, compile one dataframe for all jobs to plot.', usage='python runGREGOR.py ./index_SnpFiles/ ./index_BedFiles -t 4:00:00 -s 6 -name chromStates -ld 0.8,0.9')
    parser.add_argument('snpFileDirectory', type=str, help="""The directory with SNP files- each ending in .txt GREGOR job would be submitted for each SNP file. Should be in chr:pos or rsID format; IMPORTANT: filename should not contain '_'""")
    parser.add_argument('bedFileDirectory', type=str, help="""The directory with (unzipped) bed files. Will assume GREGOR be run using ALL files in specified directory. If GRGEOR need not be run using all files, the file bedFileIndes_xyz.txt can be manually changed before submitting drmr job.""")
    ## Provide path to sample .conf file to the 'default' variable here
    parser.add_argument('-f','--sampleConfFile', nargs='?', default='/home/arushiv/toolScripts/gregor_sampleConfFile.conf',
                        help="""Provide path to a sample .conf file where parameters will be edited as supplied to this script. (default: /lab/arushiv/toolScripts/gregor_sampleConfFile.conf)""")
    parser.add_argument('-n','--bedFilesPerJob', nargs='?', default='all', help="""Divide total bedfiles into separate GREGOR jobs. CAUTION: Only use if TOO many bed files make GREGOR job too long to run. For example, if enrichment has to be run for 3000 bed files and supplied -n is 300, this script will divide bed files into 10 jobs with 300 bedfiles each. (default: GREGOR runs all bedfiles in one job for each SNP file)""")
    parser.add_argument('-ld', '--linkageDisequilibrium', nargs='?', default='0.8', help="""minimum LD r2 for proxy SNPs. Multiple LD values to be supplied comma (,) separated. (default: 0.8)""")
    parser.add_argument('-t','--walltime', nargs='?', type=str, default='12:00:00', help="""Walltime for each GREGOR job. (default: 12:00:00 [hh:mm:ss])""")
    parser.add_argument('-s','--cores', nargs='?', type=str, default='6', help="""Number of cores for each GREGOR job. (default: 6)""")
    parser.add_argument('-d','--workingDirectory', type=str, nargs='?', default='.', help="""Directory where output directories will be created. (default: current directory)""")
    parser.add_argument('--now', action='store_const', const='now', default='wait', help="""Submit drmr job file now.""")
    parser.add_argument('-name','--nameGregorRun', type=str, nargs='?', default="", help="""Name string to be appended to conf files and output folders indicating specific gregor run. (default: "")""")
    args = parser.parse_args()

    snpFileDirectory=args.snpFileDirectory
    bedFileDirectory=args.bedFileDirectory
    sampleConfFile=args.sampleConfFile
    # bedFilesPerJob=int(args.bedFilesPerJob)
    if args.bedFilesPerJob != "all":
        bedFilesPerJob=int(args.bedFilesPerJob)
    else:
        bedFilesPerJob="all"
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


