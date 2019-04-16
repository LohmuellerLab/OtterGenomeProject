# -*- coding: utf-8 -*-
"""
Created on Mon Feb 11 15:21:42 2019

@author: annabelbeichman
"""

### ABBA-BABA for otters
# input: vcf file (remove ./. genotypes) with all three otters, mapped to outgroup
# tree: list of otters in order of outgroup > ingroup
# so giant otter, northern or southern sea otter, northern or southern sea otter
# if was apes could do gorilla, chimp, human for instance
# reference genome needs to be outgroup to other 3 populations
# and reference allele is assumed to be ancestral (this is flawed but a proxy)


import os
import sys
import argparse
import csv
import gzip
import random


parser = argparse.ArgumentParser(description='')
parser.add_argument("--vcf",required=True,help="path to vcf file that has been merged and fully filtered (no ./. genotypes; all quality filtering done). Exclude sites that are 0/0 across all individuals to increase speed. (but must have been part of initial merge so that they aren't counted as missing data in one of the individuals)")
parser.add_argument("--treeOrder",required=True,help="list of individuals in tree (excluding the reference, which must be the outgroup) in order of outgroup -> ingroup: where A is the reference genome and is unlisted. So giant otter, northern/southern sea otter, southern/northern sea otter; if was apes could do chimp, bonobo, human if mapped to a gorilla outgroup, for instance.")
parser.add_argument("--outdir",required=True,help="path to output directory")
parser.add_argument("--outPREFIX",required=False,help="output file prefix; for instance, SGE_TASK_ID (optional)",default="")

args = parser.parse_args()
vcfFile=args.vcf
treeOrder=list(args.treeOrder.split(","))
outfilepath1=args.outdir+"/ABBA-BABA."+str(args.outPREFIX)+".txt"
outfilepath2=args.outdir+"/ABBA-BABA.DStatistic."+str(args.outPREFIX)+".txt"


# So: 
# get sample names
def getSamplesFromVCF(filepath):
    inVCF = gzip.open(filepath, 'r')
    
    samples=[]
    for line in inVCF:
        if line.startswith('##'):
            pass
        else:
            for i in line.split()[9:]: samples.append(i)
            break
    inVCF.close()
    return samples
samples  = getSamplesFromVCF(vcfFile) # these are what are in the VCF file

# reset vcf:

# convert genotypes to single numbers (0 or 1); heterozygotes get randomly assigned
# with 50-50 odds
def convertGT(gt):
    if gt == "./.":
        print("Something is wrong! There are ./. genotypes present!")
        #break
    elif gt == "0|0" or gt=="0/0": 
        convGT=0 # homozygous reference
        return convGT
    elif gt=="0/1" or gt == "1/0" or gt=="1|0" or gt=="0|1":
        convGT=random.randint(0,1)
        return convGT
    elif gt=="1/1" or gt=="1|1":
        convGT=1
        return convGT
    else:
        print("A weird genotype has emerged! Check your file! "+str(gt))

# 
################# set up populations ###########
# get sample IDs from vcf and make correspond to tree order ((A=REF)B,C,D)

# set Identity of B, C, and D (in order of tree input)
popA="REFERENCE" # reference will always be 0 
popB=treeOrder[0] # first in the list given at input
popC=treeOrder[1] # second in list given at input 
popD=treeOrder[2] # third in list given at input


################## open outfile ###################
# header
outfile1=open(outfilepath1,"w")
outfile1.write("scaffold\tpos\tabba\tbaba\tA:REFERENCE\tB:"+popB+"\tC:"+popC+"\tD:"+popD+"\n")
# set up empty abbaCount and babaCount (overall Counts) for D-test
abbaCount=0
babaCount=0
################# go through vcf ##################
VCF = gzip.open(vcfFile, 'r')
for line in VCF:
    if not line.startswith("#"):
        abba=0 
        baba=0
        line = line.rstrip("\n")
        line = line.split("\t")
        myscaff=line[0] #scaffold
        mypos=line[1] # position
        # need to offset genotype fields to match with sample names
        # I do that somewhere
        mygenoinfo=line[9:] # all genotypes for all individuals
        allCalls=[i.split(":")[0] for i in mygenoinfo] # genotype calls
        callDict = dict(zip(samples,allCalls)) # zip with names
        # get converted genotypes for each population
        # input is the actual genotype (e.g. 0/0)
        # and the function converts to a single number
        # 0/0 --> 0; 1/1 --> 1 and 0/1 --> 50/50 assignment to 0 and 1
        conv_gtA=0 # ref is always ancestral
        conv_gtB=convertGT(callDict[popB])
        conv_gtC=convertGT(callDict[popC])
        conv_gtD=convertGT(callDict[popD])
        # now need to check my ABBA-BABAs
        # since reference is *always* 0 (ancestral)
        # then the only two cases we need to check are 0 1 1 0 and 0 1 0 1
        # because there will never be 1 0 0 1 or 1 0 1 0 (because A is never 1)
        if conv_gtB==1 and conv_gtC==1 and conv_gtD==0: #(0) 1 1 0 
            abba=1
            abbaCount+=1
        if conv_gtB==1 and conv_gtC==0 and conv_gtD==1: # (0) 1 0 1
            baba=1
            babaCount+=1
        # skip sites that aren't abba or baba 
        if abba==0 and baba==0:
            continue
        else:
            outfile1.write(myscaff+"\t"+str(mypos)+"\t"+str(abba)+"\t"+str(baba)+"\t0\t"+str(conv_gtB)+"\t"+str(conv_gtC)+"\t"+str(conv_gtD)+"\n")

outfile1.close()
VCF.close()
# Doing D-calc in R 
