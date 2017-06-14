import subprocess as sp
import os
import argparse
import glob
import sys
import shutil

def printWait(workflowFile):
    workflowFile.write("\n# drmr:wait\n")

def printInfo(workflowFile, info):
    workflowFile.write("\n###\n %s \n###\n" % (info))

def printResources(workflowFile, nodes, cores, memory, time):
    workflowFile.write("\n# drmr:job nodes=%s processors=%s processor_memory=%s time_limit=%s\n" %(nodes, cores, memory, time))

def printCommand(workflowFile, cmd):
    workflowFile.write('\n')
    workflowFile.write(repr(cmd).replace('"', ''))
    workflowFile.write('\n')
    
def getAllFootprintOverlaps(workflowFile, inputfile, footprintPathFile, outputfile):
    cmd = "for i in `cat %s`; do b=`basename $i .Footprints.bed | sed -e 's:Islet-union.::g' | sed -e 's/__/::/g'`; zcat /lab/work/albanus/isletEndoC_ver2/posteriors_all/bound/abcu196/ABCU196.${b}.bound.bed.gz /lab/work/albanus/isletEndoC_ver2/posteriors_all/bound/abe1388/HP14149-01.${b}.bound.bed.gz | less | awk '{print $1,$2,$3,$4,$6}' OFS='\t' | sortBed -i - | uniq | intersectBed -a %s -b - -wa -wb; done | awk '{print $1,$2,$3,$7,$8,$9,$10,$11,$5,$6}' OFS='\t' | sort | uniq > %s" % (footprintPathFile, inputfile, outputfile)     # Syntax to submit drmr job

    printCommand(workflowFile, cmd)
    
def getSNPMotifOverlap(workflowFile, outputfile):
    cmd = "python /home/arushiv/snpMotifOverlap.py %s > temp; mv temp %s" % (outputfile, outputfile)
    workflowFile.write(cmd)

    
def getAlleleFrequencies(workflowFile, outputfile, intermediateFile):
    cmd = "less %s | awk '{print $1,$2,$3}' OFS='\t' | sed -e 's:chr::g' | sort | uniq > %s; vcftools --gzvcf /lab/data/eqtl/2017_02_16_insPIRE-eQTL-vcfFiles/Islets_20Sep2016_eQTL_clean_NoMissing.vcf.gz --bed %s --freq  --out %s" %(outputfile, intermediateFile, intermediateFile, workflowName)

    printCommand(workflowFile, cmd)

def analyzeDirectionEffect(workflowFile, outputfile, freqFile, workflowOutputFile, infoContent):
    cmd = "python analyzeDirectionEffect_eqtlGwas.py %s %s %s -ic %s" %(outputfile, freqFile, workflowOutputFile, infoContent)
    workflowFile.write(cmd)
    
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Workflow for generating drmr file for getting all footprint overlap information with eQTL', usage='python workflow.py ')
    parser.add_argument('-d','--workingDirectory', type=str, nargs='?', default='.', help="""Directory where output directories will be created. (default: current directory)""")

    parser.add_argument('inputfile', type=str, help="""SNP bed file with chrom, start, end, indexSNPPos. No header""")
    parser.add_argument('-p','--footprintPathFile', type=str, default='/lab/work/arushiv/chromatin/insPIRE/new_conditional_datasets/binByBeta/enrichment_footprints/gregor_fullsetAndOverlappingBetaBinsAfterPruning/exploreEnrichedFootprints/exploreColocalizationFootprints/allTfPaths1.txt',
                        help="""Provide paths to footprint files""")

    parser.add_argument('-ic','--infoContent', type=float, default=0.7, help="""Information content filter for SNP overlap with footprint motif. Default = 0.7""")
    # parser.add_argument('-t','--walltime', nargs='?', type=str, default='12:00:00', help="""Walltime for each GREGOR job. (default: 12:00:00 [hh:mm:ss])""")
    # parser.add_argument('-s','--cores', nargs='?', type=str, default='6', help="""Number of cores for each GREGOR job. (default: 6)""")
    parser.add_argument('-name','--workflowName', type=str, nargs='?', default="workflow", help="""Name string for current workflow. Will be appended to files generated (default: "workflow")""")
    parser.add_argument('--now', action='store_const', const='now', default='wait', help="""Submit drmr job file now.""")
        
    args = parser.parse_args()

    workingDirectory = args.workingDirectory
    inputfile = args.inputfile
    footprintPathFile = args.footprintPathFile
    workflowName = args.workflowName
    workflowFile = os.path.join(workflowName+'.dr')
    infoContent = args.infoContent

    with open(workflowFile, 'w') as f:

        
        printInfo(f, "# Footprints that overlap SNPs in the input file")
        printResources(f, 1, 1, 4000, "2:00:00")
        outputfile = os.path.join(workingDirectory, workflowName + '.FootprintOverlaps.bed')
        getAllFootprintOverlaps(f, inputfile, footprintPathFile, outputfile)
        printWait(f)

        
        printInfo(f, "# Get SNP motif overlap information" )
        printResources(f, 1, 1, 4000, "1:00:00")
        getSNPMotifOverlap(f, outputfile)
        printWait(f)

        printInfo(f, "# Get allele frequencies for SNPs")
        printResources(f, 1, 1, 8000, "1:00:00")
        intermediateFile = os.path.join(workingDirectory, workflowName + ".3cols.bed")
        getAlleleFrequencies(f, outputfile, intermediateFile)
        printWait(f)

        printInfo(f, "# Get Alleles at positions of overlap")
        printResources(f, 1, 1, 4000, "1:00:00")
        freqFile = os.path.join(workingDirectory, workflowName + '.frq')
        workflowOutputFile = os.path.join(workingDirectory, workflowName + '.finalOutput.dat')
        analyzeDirectionEffect(f, outputfile, freqFile, workflowOutputFile, infoContent)


    if args.now != "wait":
        sp.call("drmr %s" % workflowFile, shell=True)
