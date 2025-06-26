import vcf
import sys
import argparse
import numpy as np
import pandas as pd

#Parse arguments
parser = argparse.ArgumentParser(description='Information about number of sliding windows and step size')
parser.add_argument('-vcf_file', dest = 'vcf_file', action='store', nargs = 1, type = str, help = 'path to vcf')
parser.add_argument('-outFile', dest = 'outFile', action='store', nargs = 1, type = str, help = 'path to output file (suffix will be added)')
parser.add_argument('-samples', dest = 'samples', action='store', nargs = 1, type = int, help = 'number  of chromosomes sampled')
parser.add_argument('-ro', dest = 'ro', action='store', nargs = 1, type = float, help = 'recombination rate')
parser.add_argument('-sample_type', dest = 'sample_type', action='store', nargs = 1, type = str, help = 'whether sample is chimeric or not')

#read input parameters
args = parser.parse_args()
vcf_file = args.vcf_file[0]
out_file = args.outFile[0]
samples = args.samples[0]
ro = args.ro[0]
sample_type = args.sample_type[0]



#Get chimeric site counts from vcf file
def get_chimeric_aff(vcf_file, samples, ro):
    vcf_reader = vcf.Reader(filename=vcf_file)
    d_af = {}
    for i,record in enumerate(vcf_reader):
        #create list of genotypes
        lst = (''.join([x['GT'] for x in record.samples]).replace('|',''))
        #isolate non-identeical twins
        twins = [lst[x:x+4] for x in range(0, len(lst), 4)]
        #Loop through twins, identifying number of derived sites
        c=0
        for t in twins:
            #Check first and third sites, and then second and fourth (ie match chromosomes)
            if((t[0]=='1') | (t[2]=='1')):
                c+=1
            if((t[1]=='1') | (t[3]=='1')):
                c+=1
        d_af[record.POS] = c

    df = pd.DataFrame([d_af.keys(),d_af.values()]).T
    df.columns = ['physPos','x']
    df = df[df['x']>0]
    df['n'] = samples

    #Add ro info
    df['genPos'] = df.physPos * ro
    df = df[['physPos', 'genPos', 'x', 'n']]
    df = df.drop_duplicates()

    aff = df
    aff['physPos'] = [int(x) for x in aff.physPos]
    aff['x'] = [int(x) for x in aff['x']]
    aff['n'] = [int(x) for x in aff['n']]
    return(aff)


#Get chimeric site counts from vcf file
def get_WF_aff(vcf_file, samples, ro):
    vcf_reader = vcf.Reader(filename=vcf_file)
    d_af = {}
    for i,record in enumerate(vcf_reader):
        #create list of genotypes
        lst = (''.join([x['GT'] for x in record.samples]).replace('|',''))
        c = lst.count('1')
        d_af[record.POS] = c

    df = pd.DataFrame([d_af.keys(),d_af.values()]).T
    df.columns = ['physPos','x']
    df = df[df['x']>0]
    df['n'] = samples

    #Add ro info
    df['genPos'] = df.physPos * ro
    df = df[['physPos', 'genPos', 'x', 'n']]
    df = df.drop_duplicates()

    aff = df
    aff['physPos'] = [int(x) for x in aff.physPos]
    aff['x'] = [int(x) for x in aff['x']]
    aff['n'] = [int(x) for x in aff['n']]
    return(aff)


if(sample_type=='chimeric'):
    aff = get_chimeric_aff(vcf_file, samples, ro)
else:
    aff = get_WF_aff(vcf_file, samples, ro)

aff.to_csv(out_file + '.aff', sep='\t', header=True, index=False)

