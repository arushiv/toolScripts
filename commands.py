
import subprocess as sp
import os
import argparse
import glob
import sys
import shutil
import pandas
import csv
import errno

def printWait(workflowFile):
    workflowFile.write("\n# drmr:wait\n")

def printInfo(workflowFile, info):
    workflowFile.write("\n###\n # %s \n###\n" % (info))

def printResources(workflowFile, nodes, cores, memory, time):
    workflowFile.write("\n# drmr:job nodes=%s processors=%s processor_memory=%s time_limit=%s\n" %(nodes, cores, memory, time))

def withIonIce(x):
    return "ionice -c2 -n7 " + x

def newmkdir(dirname):
    try:
        os.makedirs(dirname)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(dirname):
            pass
        else:
            raise

def printBash(x):
    return x.replace('\t', '\\t')

    
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Workflow for generating drmr file for ', usage='python commands.py ')
    parser.add_argument('-d','--workingDirectory', type=str, nargs='?', default='.', help="""Directory where output directories will be created. (default: current directory)""")
    parser.add_argument('-name','--workflowName', type=str, nargs='?', default="workflow", help="""Name string for current workflow. Will be appended to files generated (default: "workflow")""")
    parser.add_argument('--now', action='store_const', const='now', default='wait', help="""Submit drmr job file now.""")
    # parser.add_argument('-f', '--footprintDirPath', type=str, default='/lab/work/arushiv/chromatin/crossSpeciesAnalysis/mappedToHg19Footprints/hmr05Subsetted', help="""Path to directory containing footprints""")
    # parser.add_argument('-m', '--motifFilePath', type=str, default='/lab/work/arushiv/chromatin/crossSpeciesAnalysis/mappedToHg19Footprints/allCellNonEmptyFootprints.hmr05subsetted.dat', help="""Path to directory containing footprints""")

    args = parser.parse_args()

    workingDirectory = args.workingDirectory
    intermediateFileDir = os.path.join(workingDirectory, "intermediateFiles")
    newmkdir(intermediateFileDir)
    workflowName = args.workflowName
    workflowFile = os.path.join(workflowName+'.dr')


    with open(workflowFile, 'w') as f:

        
        printInfo(f, "# ")
        printResources(f, 1, 1, 4000, "30:00")

        printWait(f)

        printInfo(f, "# ")
        printResources(f, 1, 1, 4000, "10:00")
        outputfile = os.path.join(intermediateFileDir, "gatResults.dat")

        printWait()
        
        printInfo(f, "#")
        printResources(f, 1, 1, 4000, "10:00")
        figfile = os.path.join(intermediateFileDir, "fig.2wayGat.pdf")
        

    if args.now != "wait":
        sp.call("drmr %s" % workflowFile, shell=True)
