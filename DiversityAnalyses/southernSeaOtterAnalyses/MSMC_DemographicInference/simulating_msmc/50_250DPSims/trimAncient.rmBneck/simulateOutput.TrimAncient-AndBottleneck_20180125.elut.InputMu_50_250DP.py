# -*- coding: utf-8 -*-
"""
Created on Tue Oct 31 13:40:10 2017

@author: annabelbeichman
"""
# This is to set up a MACS script to simulate the full demographic history of the sea otter or giant otter 
# want to make generic
# input
# spp
# mu
# Len

import pandas as pd

import argparse
parser = argparse.ArgumentParser(description='Make a macs simulation from msmc')
parser.add_argument("--mu",required=True,help="supply mutation rate in mutation/bp/gen")
args = parser.parse_args()
mu= float(args.mu)


trimIndex=33 # 0-based index of where you want to trim
r = 1e-08
Len = 30000000
num=60
blocksPerGroup=10
groups=num/blocksPerGroup
input=pd.read_table("/Users/annabelbeichman/Documents/UCLA/Otters/RawReadsToGenotypes_Dec2016/mappingGidgetReadsToFerret/MSMC_DemographicInference_WholeGenome/msmc/output_20180209_50_250DPFilter/sea_otter.msmc.out.final.txt")
input.head()
diploids=tuple((1/x)/(2*mu) for x in input.lambda_00)
diploids_trimancient = diploids[0:trimIndex]
#### Want to try eliminating the most recent bottleneck ; doesn't affect Na or anything

# So want that 526 to just stay the same as 1314 for that time point:
# need to make it a list to change a value:
diploids_trimancient_list = list(diploids_trimancient)
diploids_trimancient_list[1] = diploids_trimancient_list[2] # getting rid of the bottleneck and just keeping the population size constant through that time point
# now just reset it to be a tuple again with the bottleneck removed:
diploids_trimancient = tuple(diploids_trimancient_list)
### NOTE THAT THIS DOESN'T CHANGE THE TIME POINT STRUCTURE AT ALL. IT JUST REMOVES THE DIP
### can then resume the rest of the script in the same way
Na = diploids_trimancient[-1] # ancestral size 
# scale it relative to oldest time in inference (could also do most recent size, just be consistent)
diploids_trimancient_Na = diploids_trimancient/Na # this is now relative to ancient size, so is ready for macs
# now the [1] and [2] entries are the same: [ 0.19555707, [this used to be the bottleneck] 0.10518316,  0.10518316

times_gen = tuple(x/mu for x in input.left_time_boundary) # this gives time in generations (if you wanted years would multiply by X yrs/gen)
times_gen_trimancient= times_gen[0:trimIndex]
times_gen_trimancient_4Na = times_gen_trimancient / (4*Na)
# if first time is -0, make it just 0.
times_gen_trimancient_4Na[0]=0 # set first time to zero if its -0. 

ss=2 # two haplotypes (one genome)
theta = 4*Na*mu
rho = 4 * Na*r

# Note that theta will be the same across different values of Mu for the same
# models, because mutation rate changes scaling of MSMC and so alters Na, but then
# theta = 4Namu = 4 (1/(lamba*2mu))*mu = 2/lamba [is this right?]
############################## WRITE SCRIPT ################
print("#!/bin/bash")
print("#$ -cwd")
print("#$ -l h_rt=5:00:00,h_data=28G,highp")
print("#$ -N elutMSMCTrim")
print("#$ -m abe")
print("#$ -M ab08028")



print("rundate=`date +%Y%m%d`")


print("model=elut.msmc.trimAncient.rmBneck."+str(mu))
print("mkdir -p ${model}")
print("for j in {1.."+str(groups)+"}") 
print("do")
print("mkdir -p ${model}/group_$j.${model}")

print("cd ${model}/group_$j.${model}")
print("cp ../../../macs ./")
print("cp ../../../msformatter ./")
print("for i in {1.."+str(blocksPerGroup)+"}")
print("do")

print("# sea otter msmc model -- full msmc")
print("mu="+str(mu))
print("r="+str(r))
print("Na="+str(Na))
print("rho=" +str(rho))
print("theta="+str(theta))
print("date=`date +%Y%m%d`")
print("SEED=$((date+$RANDOM+((j-1)*"+str(blocksPerGroup)+")+i))") #
print("# this is a new addition! need to have a different random seed for each simulation; if they start within a second of each other, they will have the same seed. not an issue for big simulations of 30Mb because those are slow, but 100kb can start within a second of each other!")
print("./macs " +str(ss) +" "+str(Len)+" -t "+str(theta)+" -r "+str(rho))+" -s $SEED",
for x, y in zip(times_gen_trimancient_4Na,diploids_trimancient_Na):
    print("-eN " + str(x)+" "+str(y)),
print(" > group_${j}_block_${i}.${model}.macsFormat.OutputFile.${rundate}.txt"),

print("")
print("./msformatter < group_${j}_block_${i}.${model}.macsFormat.OutputFile.${rundate}.txt > group_${j}_block_${i}.${model}.msFormat.OutputFile.${rundate}.txt")
###################################################
print("done")
print("cd ../../")
print("done")
