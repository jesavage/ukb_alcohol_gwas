#!/usr/bin/env python

#
# This example simulates a case-control sample using a GxE model.
#

import sys, os, random, math, logging

import simuOpt
simuOpt.setOptions(gui=False, alleleType='binary', optimized=True, numThreads=30)
from simuPOP import *

from simuPOP.utils import importPopulation, export    
from simuPOP.utils import saveCSV

import loadHapMap3, selectMarkers, simuGWAS

def downloadData(chroms, logger):
    '''
    Download and create populations from the third phase of the HapMap3 data.
    This equivalent to command

    > loadHapMap3.py --chroms='[2, 5,10]' --dest=HapMap
    '''
    if not os.path.isdir('HapMap'):
        os.mkdir('HapMap')
    for chrom in chroms:
        for popName in loadHapMap3.HapMap3_pops:
            filename = 'HapMap/HapMap3_%s_chr%d.pop' % (popName, chrom)
            if not os.path.isfile(filename):
                pop = loadHapMap3.loadHapMapPop(chrom, popName, logger)
                pop.save(filename)

def getInitPop(logger):
    '''
    Select 2000 markers on a random regions on chromosomes 2, 5 and 10, using markers from the Illumina 1M chipset
    
    From command line, you could prepare a marker list file from illumina annotation file
         > cut -d, -f2 HumanHap550v3_A.csv  > HumanHap550v3_A.lst
    and then select markers
         > selectMarkers.py --markerList='HumanHap550v3_A.lst' --chroms='[2, 5,10]' \
            --numMarkers='[2000,2000,2000]' --startPos='[20000000, 20000000,40000000]' \
            --filename=ex2_init.pop
    '''
    if os.path.isfile('ex2_init2001.pop') and os.path.isfile('ex2_init.pop.lst'):
        if logger:
            logger.info('ex2_init2001.pop already exists. Please remove this file if you would like to regenerate an initial population.')
        return
    ann = open('hh1mv3_snptable.txt')
    names = []
    for line in ann:
        names.append(line.split(',')[0])
    if logger:
        logger.info('Select 1500 markers from chromosomes 2, 5 and 10')
    pop = selectMarkers.getHapMapMarkers(
        names=names,
        HapMap_dir='HapMap',
        chroms=[2, 5, 10],
        HapMap_pops=selectMarkers.HapMap3_pops,
        startPos=[23000000, 23000000, 44000000],
        numMarkers=[500, 500, 500],
        mergeSubPops=False,
        logger=logger)
    if logger:
        logger.info('Saving initial population to ex2_init.pop')
    pop.save('ex2_init2001.pop')
    if logger:
        logger.info('Saving marker information to ex2_init.pop.lst')
    selectMarkers.saveMarkerList(pop, 'ex2_init2001.pop.lst', logger)


## Three SNPs chosen from the middle of each region with roughly equivalent and not too small MAFs
## Assume no selection pressure (fitness values =1 for each genotype)
## Assume current allele frequencies are similar to ancestral
## Environment interaction modeled with a) strong, b) weak, c) no main SNP effect
# rs2044766	2	23179561	C	T	0.3112
# rs1521011	5	23658249	C	T	0.3651
# rs7908745	10	45273773	A	G	0.2779


DPL=["rs2044766", "rs1521011", "rs7908745"]

def expandPop(logger):
    # This is equivalent to
    #
    #  > simuGWAS.py --initPop=ex2_init.pop --DPL='["rs2044766", "rs7718389", "rs2279433"]' \
    #  --filename=gxe_expanded.pop --curAlleleFreq='[0.05, 0.15]' --trajectory='forward'  --mlSelModel='multiplicative' \
    #  --scale=1 --optimized --gui=False --fitness='[1, 0.996, 0.994, 1, 1.001, 1.005]'
    #
    # This just to make this result reproducible.
    getRNG().set(seed=1355)
    #
    filename = 'gxe_expanded2001.pop'
    if os.path.isfile(filename):
        if logger:
            logger.info('%s already exists. Please remove this file if you would like to regenerate an expanded population.' % filename)
        return
    else:
        if logger:
            logger.info('Simulating an expanded population from ex2_init.pop...')
    pop = loadPopulation('ex2_init2001.pop')
    pop.setInfoFields(['ind_id','father_id','mother_id'])
    pars = simuOpt.Params(simuGWAS.options, DPL=DPL,
            curAlleleFreq=[0.31,.37,.28], trajectory='forward',  
        trajPlot='gxe_traj.pdf', mlSelModel='additive',
        scale=1, fitness=[1,1,1 ])
    pop = simuGWAS.simuGWAS(pars, pop, logger=logger)
    if logger:
        logger.info('Saving expanded population to gxe_expanded.pop')
    pop.save('gxe_expanded2001.pop')        


selectedCase = 0
selectedControl = 0
discardedInds = 0

alpha = -2
beta1 = 0.3
beta2 = 0.1
beta3 = 0
betaE = 0.1 
gamma1 = 0.3
gamma2 = 0.3
gamma3 = 0.3

env = []
phen = []

random.seed(20001)  # to keep result reproducible.
getRNG().set(seed=2001)

def _selectInds(off, param):
    'Determine if the offspring can be kept.'
    e = random.randint(0, 1)
    g1 = off.allele(param[0], 0) + off.allele(param[0], 1)
    g2 = off.allele(param[1], 0) + off.allele(param[1], 1)
    g3 = off.allele(param[2], 0) + off.allele(param[2], 1)
    logit = alpha + beta1*g1 + beta2*g2 + beta3*g3 + betaE*e + gamma1*e*g1 + gamma2*e*g2 + gamma3*e*g3
    affected = random.random() < (1 / (1. + math.exp(-logit)))
    global selectedCase, selectedControl, discardedInds, env, phen
    if affected:
        if selectedCase < 50000:
            off.setAffected(True)
            #off.setIndInfo(e,'env')
            selectedCase += 1
            env.append(e)
            phen.append(affected)
            return True
    else:
        if selectedControl < 50000:
            selectedControl += 1
            off.setAffected(False)
            #off.setIndInfo(e,'env')
            env.append(e)
            phen.append(affected)
            return True
    discardedInds += 1
    return False

def generateSample(numCase=50000, numCtrl=50000, logger=None):
    if logger:
        logger.info('Generating %d cases and %d controls...' % (numCase, numCtrl))
    pop = loadPopulation('gxe_expanded2001.pop')
    loci = pop.lociByNames(DPL)
    pop.evolve(
        matingScheme=RandomMating(
            ops=[
                PedigreeTagger(),
                IdTagger(),
                MendelianGenoTransmitter(),
                # an individual will be discarded if _selectInds returns False
                PyOperator(func=_selectInds, param=loci)
            ], 
            subPopSize=100000
        ),
        gen = 1
    )
    pop.save('gxe_sample2001.pop')
    pop.asPedigree()
    export(pop, format='PED', output='gxe_sample2001.ped')
    if logger:
        logger.info('Disease prevalence: %.4f' % (50000. / (100000 + discardedInds)))


if __name__ == '__main__':
    # You can change logging level to DEBUG to get more information
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('example2')
#    downloadData([2, 5, 10], logger)
    getInitPop(logger)
    expandPop(logger)
    generateSample(logger=logger)
    #save phenotype variable per individual
    with open(r'gxe_sample2001.pheno', 'w') as fp:
        fp.write("\n".join(str(item) for item in phen))
    #save environmental variable per individual
    with open(r'gxe_sample2001.env', 'w') as fp:
        fp.write("\n".join(str(item) for item in env))

