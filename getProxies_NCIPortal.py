#!/usr/bin/env python3


import argparse
import errno
import glob
import os
import shutil
import subprocess as sp
import sys
import time

from bs4 import BeautifulSoup, SoupStrainer
import pandas
import requests
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as expect


def newmkdir(dirname):
    try:
        os.makedirs(dirname)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(dirname):
            pass
        else:
            raise

def getPopulationString(population):
    outlist = []
    for name in population:
        if name not in ['AFR', 'AMR', 'EAS', 'EUR', 'SAS']:
            outlist.append(name)
        elif name == "EUR":
            outlist.extend(['CEU','TSI','FIN','GBR','IBS'])
        elif name == "AFR":
            outlist.extend(['YRI','LWK','GWD','MSL','ESN','ASW','ACB'])
        elif name == "EAS":
            outlist.extend('CHB','JPT','CHS','CDX','KHV')
        elif name == "AMR":
            outlist.extend(['MXL','PUR','CLM','PEL'])
        elif name == "SAS":
            outlist.extend(['GIH' ,'PJL', 'BEB', 'STU', 'ITU'])

    return '%2B'.join(outlist)


def getProxies(outdir, populationString, variant):
    print('Fetching proxies for variant {}'.format(variant['snp']))

    link = "https://analysistools.nci.nih.gov/LDlink/?var=%s&pop=%s&r2_d=r2&tab=ldproxy" % (variant['snp'], populationString)

    driver = webdriver.PhantomJS()
    driver.set_window_size(2560, 1440)
    driver.get(link)

    # wait for results
    resultlink = WebDriverWait(driver, 600).until(
        expect.visibility_of_element_located((By.ID, "ldproxy-results-link"))
    )

    if not resultlink:
        print('Could not get proxies from {}'.format(link), file=sys.stderr)
        sys.exit(1)

    result_url = resultlink.get_attribute('href')
    r = requests.get(result_url)
    if r.status_code != 200:
        print('Could not read list of proxies from {}: {}.'.format(result_url, r.status_code), file=sys.stderr)
        sys.exit(1)

    with open(os.path.join(outdir, '{}.proxies.txt'.format(variant['snp'])), 'wb') as out:
        out.write(r.content)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Get proxy variants from NCI website https://analysistools.nci.nih.gov ', usage='python getProxies_NCIPortal.py <variantfile> -p EUR -ld 0.8 -o <outputfilename> ')
    parser.add_argument('variantfile', type=argparse.FileType('r'), help="""file with rsIDs to search proxies of.""")
    parser.add_argument('-p', '--population', action='append', default=['EUR'], help="""Population to search proxies in. Options are: \n
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
    parser.add_argument('-r', '--resultFileDir', default='results_proxySearch', help="""Directory that will contain all proxy search results.""")
    # parser.add_argument('outputfile', help="""name of the outputfile""")
    parser.add_argument('-ld', '--ldThreshold', type=float, default=0.8, help="""Fetch proxies with this r2 or higher and paste into the output file""")

    args = parser.parse_args()
    newmkdir(args.resultFileDir)

    population = args.population
    populationString =  getPopulationString(population)

    d = pandas.read_csv(args.variantfile, sep='\t', usecols = ['snp'])

    d.apply(lambda x: getProxies(args.resultFileDir, populationString, x), axis=1)
