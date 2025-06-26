import pandas as pd
import numpy as np
import random
import argparse
import sys

#parsing user given constants
parser = argparse.ArgumentParser(description='Information about length of region and sample size')
parser.add_argument('-pedigrees', dest = 'pedigrees', action='store', nargs = 1, type = str, help = 'path to pedigree file')
parser.add_argument('-outFile', dest = 'outFile', action='store', nargs = 1, type = str, help = 'path to store output files')
parser.add_argument('-samples', dest = 'samples', action='store', nargs = 1, type = int, help = 'no. of samples required')

args = parser.parse_args()
samples = args.samples[0]
pedigrees = args.pedigrees[0]
outFile = args.outFile[0]

df = pd.read_csv(pedigrees, sep='\t', names=['id', 'pedigreeID', 'parent1', 'parent2'])
#Check if sample is chimeric by identifying duplicates (ie same parents non-identical twins)
m1 = df.duplicated(['parent1','parent2'], keep=False)
df['chimerism'] = np.select([m1],[1], default=0)

#Get non-chimeric samples and output to file
tdf = df[df.chimerism==0]
p1 = list(tdf['parent1'])
p2 = list(tdf['parent2'])

#lst = random.sample(range(0, len(p1)-1), samples)
#l1 = [p1[i] for i in lst]
#l2 = [p2[i] for i in lst]
#u = pd.DataFrame({"parent1":l1,"parent2":l2})
#rdf = df.merge(u)
#rdf[['id']].to_csv(outFile + "_non-chimeric.txt", header=False, index=False)

#Repeat for chimeric samples
tdf = df[df.chimerism==1].drop_duplicates(['parent1', 'parent2'])
p1 = list(tdf['parent1'])
p2 = list(tdf['parent2'])
lst = random.sample(range(0, len(p1)-1), samples)
l1 = [p1[i] for i in lst]
l2 = [p2[i] for i in lst]
u = pd.DataFrame({"parent1":l1,"parent2":l2})
rdf = df.merge(u)
rdf[['id']].to_csv(outFile + "_chimeric.txt", header=False, index=False)