#!/usr/bin/env python3.5

import subprocess as sp
import os
import glob
import sys
import pandas
from snakemake.utils import R

DIRECTORIES = {
    'data' : "data",
    'intermediateFiles' : "intermediateFiles",
    'scripts' : "scripts",
    'figures' : "figures"  
}    


rule final:
    input:
        os.path.join(DIRECTORIES['figures'], "fig.deviationFromExpectation.pdf")

rule makeDirectories:
    output:
            DIRECTORIES['data'],
            DIRECTORIES['intermediateFiles'],
            DIRECTORIES['scripts'],
            DIRECTORIES['figures']
    shell:
        r"""mkdir -p {output}"""

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
