
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

def printGatCommands(workflowFile, motifFilePath, footprintDirPath, gatOutput):
    def iterateCells(x):
    # MIN6.ARNTL_1.hg19Mapped_150ext_hmr05subsetted.Footprints.bed
        cells = [['Islet-union','hmr05subsetted.Footprints.bed'],['endoC','hmr05subsetted.Footprints.bed'],['MIN6', 'hg19Mapped_150ext_hmr05subsetted.Footprints.bed'],['INS1','hg19Mapped_150ext_hmr05subsetted.Footprints.bed'],['RI','hg19Mapped_150ext_hmr05subsetted.Footprints.bed']]
        commands = ""
        for cell1, cell2 in zip(cells[0:4], cells[1:5]):
            file1 = os.path.join(footprintDirPath, "%s.%s.%s" % (cell1[0], x, cell1[1]))
            file2 = os.path.join(footprintDirPath, "%s.%s.%s" %(cell2[0], x, cell2[1]))
            outputfile = os.path.join(gatOutput, "%s.%s.%s.out" %(x, cell1[0], cell2[0]))
            cmd = "gat-run.py -a %s -s %s -w /lab/work/arushiv/chromatin/crossSpeciesAnalysis/mappedToHg19Footprints/hmr_tiletrack.0.5.concatenated.bed  | grep -v \"#\" > %s; " %(file1, file2, outputfile)
            commands = commands + withIonIce(cmd)
        return commands

    motifFile = pandas.read_csv(motifFilePath, sep='\t', usecols=['motif'])
    cmddf = motifFile.motif.map(lambda x: iterateCells(x)).to_frame()
    cmddf.to_csv(workflowFile, header=None, index=None, sep=' ', mode='a', quoting=csv.QUOTE_NONE, escapechar=' ')
    
def analyzeGATResults(workflowFile, gatOutput, outputfile):
    cmd = "python ~arushiv/toolScripts/analyze_GATResults.py -s *.out -d %s %s --split -ic motif cell1 cell2" %(gatOutput, outputfile)
    workflowFile.write(cmd)

def plot(workflowFile, outputfile, figfile):
    cmd = "Rscript plot.R %s %s" %(outputfile, figfile)
    workflowFile.write(cmd)
    
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Workflow for generating drmr file for ', usage='python commands.py ')
    parser.add_argument('-d','--workingDirectory', type=str, nargs='?', default='.', help="""Directory where output directories will be created. (default: current directory)""")
    parser.add_argument('-name','--workflowName', type=str, nargs='?', default="workflow", help="""Name string for current workflow. Will be appended to files generated (default: "workflow")""")
    parser.add_argument('--now', action='store_const', const='now', default='wait', help="""Submit drmr job file now.""")
    parser.add_argument('-f', '--footprintDirPath', type=str, default='/lab/work/arushiv/chromatin/crossSpeciesAnalysis/mappedToHg19Footprints/hmr05Subsetted', help="""Path to directory containing footprints""")
    parser.add_argument('-m', '--motifFilePath', type=str, default='/lab/work/arushiv/chromatin/crossSpeciesAnalysis/mappedToHg19Footprints/allCellNonEmptyFootprints.hmr05subsetted.dat', help="""Path to directory containing footprints""")

    args = parser.parse_args()

    workingDirectory = args.workingDirectory
    intermediateFileDir = os.path.join(workingDirectory, "intermediateFiles")
    newmkdir(intermediateFileDir)
    workflowName = args.workflowName
    workflowFile = os.path.join(workflowName+'.dr')

    footprintDirPath = args.footprintDirPath
    motifFilePath = args.motifFilePath

    gatOutput = os.path.join(workingDirectory, workflowName + "_gatOutputs")

    with open(workflowFile, 'w') as f:

        
        printInfo(f, "# Run GAT")
        printResources(f, 1, 1, 4000, "30:00")
        printGatCommands(f, motifFilePath, footprintDirPath, gatOutput)
        printWait(f)

        printInfo(f, "# Compile GAT Results")
        printResources(f, 1, 1, 4000, "10:00")
        outputfile = os.path.join(intermediateFileDir, "gatResults.dat")
        printGatCommands(f, gatOutput, outputfile)
        printWait()
        
        printInfo(f, "# Plot")
        printResources(f, 1, 1, 4000, "10:00")
        figfile = os.path.join(intermediateFileDir, "fig.2wayGat.pdf")
        printGatCommands(f, outputfile, figfile)


    if args.now != "wait":
        sp.call("drmr %s" % workflowFile, shell=True)
