
import argparse
import errno
import os
import subprocess as sp
import glob
import pandas

"""
Script to get variants in LD using 1000g phase 3 data using the NCI portal https://analysistools.nci.nih.gov/LDlink/?tab=ldproxy Can also provide a MAF threshold to filter variants.
"""
## RS_Number Coord Alleles MAF Distance Dprime R2 Correlated_Alleles RegulomeDB Function

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
            print("!!! result directory exists. CHECK COMMAND !!!")
            pass
        else:
            raise

def getOpts():
    parser = argparse.ArgumentParser(description='Script to get LD proxies from specific populations using 1000g phase 3 data using the NCI portal https://analysistools.nci.nih.gov/LDlink/?tab=ldproxy.', usage='python pruneVariantsAndFilterMaf_1000gNCIPortal.py <rsIDlist> -p EUR -r2 0.8 <outputfilename> ')
    parser.add_argument('variantfile', type=argparse.FileType('r'), help="""file with rsIDs to search proxies of. Also needs chrom information. Tab separated. Header name of the rsID and chrom column should be 'snp', 'chrom'. Pre-sort this by pre-value etc, variants further down in the list that are in high LD with top variants will be removed""")
    parser.add_argument('drmrfile', help="""name of the drmr file""")
    parser.add_argument('outputfile', help="""name of the outputfile""")
    parser.add_argument('-hv', '--variantfileHeader', nargs='+', help="""If variant file does not contain header, supply space separated list here. Header name of the rsID and chrom column should be 'snp', 'chrom'""")
    parser.add_argument('-p', '--population', default='EUR', help="""Population to search proxies in. Options are: \n
(AFR) African:[[
    (YRI) Yoruba in Ibadan, Nigera
    (LWK) Luhya in Webuye, Kenya
    (GWD) Gambian in Western Gambia
    (MSL) Mende in Sierra Leone
    (ESN) Esan in Nigera
    (ASW) Americans of African Ancestry in SW USA
    (ACB) African Carribbeans in Barbados ]]
(AMR) Ad Mixed American [[
    (MXL) Mexican Ancestry from Los Angeles, USA
    (PUR) Puerto Ricans from Puerto Rico
    (CLM) Colombians from Medellin, Colombia
    (PEL) Peruvians from Lima, Peru ]]
(EAS) East Asian [[
    (CHB) Han Chinese in Bejing, China
    (JPT) Japanese in Tokyo, Japan
    (CHS) Southern Han Chinese
    (CDX) Chinese Dai in Xishuangbanna, China
    (KHV) Kinh in Ho Chi Minh City, Vietnam ]]
(EUR) European [[
    (CEU) Utah Residents from North and West Europe
    (TSI) Toscani in Italia
    (FIN) Finnish in Finland
    (GBR) British in England and Scotland
    (IBS) Iberian population in Spain ]]
(SAS) South Asian [[
    (GIH) Gujarati Indian from Houston, Texas
    (PJL) Punjabi from Lahore, Pakistan
    (BEB) Bengali from Bangladesh
    (STU) Sri Lankan Tamil from the UK
    (ITU) Indian Telugu from the UK \n ]]
    provide sub or superpopulation code(s) WITHOUT parantheses.""" )
    parser.add_argument('-dir', '--resultFileDir', default='results_ldPrune', help="""Directory that will contain all LD prune results.""")
    parser.add_argument('-r2', '--r2Threshold', type=float, default=0.99, help="""Retain variants with this r2 or lower and paste into the output file. Default=0.99""")
    parser.add_argument('--now', action='store_const', const='now', default='wait', help="""Submit drmr job file now.""")    
    return parser

def getProxyNCI(workflowFile, d, resultFileDir, population, r2Threshold, outputfile):
    parameters = {
        'population' : population,
        'r2Threshold' : r2Threshold,
    }

    def printCommands(x):
        new = {
            'variant' : x,
            'result' : os.path.join(resultFileDir, "%s_ProxyResult.dat.gz" %(x))
        }
        parameters.update(new)
            
        cmd = """curl -k -X GET 'https://analysistools.nci.nih.gov/LDlink/LDlinkRest/ldproxy?var={variant}&pop={population}&r2_d={r2Threshold}' | gzip -c > {result} \n""".format(**parameters)
        workflowFile.write(cmd)

    d.drop_duplicates('snp')['snp'].map(lambda x: printCommands(x))

    printWait(workflowFile)
    printResources(workflowFile, 1, 1, 4000, "15:00")
    
    cmd_mergeOutputs = r"""for i in `ls {resultfiles}`; do b=`basename $i _ProxyResult.dat`; zcat $i | awk -F'\t' '{{if (( $7>={r2Threshold} )) print $1,$2,$3,$4,$5,$7"\t""'"$b"'"}}' OFS='\t'; done | sort | uniq > {output}; echo -e "proxy_rsID\tproxy_Coord\tproxy_Alleles\tproxy_MAF\tdistance\tR2\tsnp" | cat - {output} > temp; mv temp {output}""".format(resultfiles=os.path.join(resultFileDir, "*_ProxyResult.dat.gz"), r2Threshold=r2Threshold, output=outputfile)
    workflowFile.write(cmd_mergeOutputs)
                    

if __name__ == '__main__':
    parser = getOpts()
    args = parser.parse_args()
    resultFileDir = args.resultFileDir
    newmkdir(args.resultFileDir)
    variantfile = args.variantfile
    population = args.population
    r2Threshold = args.r2Threshold
    outputfile = args.outputfile
    workflowFile = args.drmrfile

    if args.variantfileHeader is not None:
        d = pandas.read_csv(variantfile, sep='\t', header=None, names=args.variantfileHeader)
    else:
        d = pandas.read_csv(variantfile, sep='\t')

    with open(workflowFile, 'w') as f:
        printInfo(f, "Get LD proxies using 1000g VCF files from NCI portal")
        printResources(f, 1, 1, 4000, "45:00")
        getProxyNCI(f, d, resultFileDir, population, r2Threshold, outputfile)

    if args.now != "wait":
        sp.call("drmr %s" % workflowFile, shell=True)
