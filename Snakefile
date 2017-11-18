import os
import pandas

DATA = {
    'summaryStats' : "/lab/data/gwas/2017_09_MAGIC/summaryStatistics/MAGIC_1KG_EA_{trait}_InvNorm_meta_2GC.gz",
    'chromStateDir' : "/lab/work/arushiv/newhuman_datasets_by_13chromatinStates/Islets.{chromState}.bed"
}

DIRECTORIES = {
    'data' : "data",
    'intermediateFiles' : "intermediateFiles",
    'scripts' : "scripts",
    'figures' : "figures"  
}    


rule final:
    input:
        os.path.join(DIRECTORIES['figures'], "fig.deviationFromExpectation.pdf")

rule linkDataFile:
    input:
        <>
    output:
        os.path.join(DIRECTORIES['data'], "<>")
    shell:
        "ln -s {input} {output}"

rule plot:
    input:
        rules.intersectAnnotations.output
    output:
        os.path.join(DIRECTORIES['figures'], "fig.deviationFromExpectation.pdf")
    script:
        os.path.join(DIRECTORIES['scripts'], "plot.R")
