import vcf
import sys
import dadi
import argparse

#SCRIPT TO OBTAIN SFS FROM VCF FILE. EXAMPLE INPUT:
#demog=(expansion eq contraction bneck)
#for d in "${demog[@]}"; do for i in $(seq 1 100); do python3 get_sfs_for_fsc.py -vcf_file ../results/demog_only/"$d"/rr_fixed_mu_fixed/rep"$i"_chimeric.vcf.gz \ 
#-outFile ../fsc/sfs_files/demog_only/"$d"/rr_fixed_mu_fixed/rep"$i"_chimeric_DAFpop0.obs -region_length 100000 -samples 20 -sample_type chimeric; done; done

#Parse arguments
parser = argparse.ArgumentParser(description='Information about number of sliding windows and step size')
parser.add_argument('-vcf_file', dest = 'vcf_file', action='store', nargs = 1, type = str, help = 'path to vcf')
parser.add_argument('-outFile', dest = 'outFile', action='store', nargs = 1, type = str, help = 'path to output file')
parser.add_argument('-samples', dest = 'samples', action='store', nargs = 1, type = int, help = 'number  of individuals for which the SFS will be made')
parser.add_argument('-sample_type', dest = 'sample_type', action='store', nargs = 1, type = str, help = 'whether sample is chimeric or not')

#read input parameters
args = parser.parse_args()
vcf_file = args.vcf_file[0]
outFile = args.outFile[0]
samples = args.samples[0]
sample_type = args.sample_type[0]


#Get chimeric site counts from vcf file
def get_chimeric_counts(vcf_file):
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
        d_af[i] = c
    
    l_af = [d_af[x] for x in d_af.keys()]
    return(l_af)


#Get non-chimeric site counts from vcf file
def get_counts(vcf_file):
    vcf_reader = vcf.Reader(filename=vcf_file)
    d_af = {}
    for i,record in enumerate(vcf_reader):
        d_af[i] = (''.join([x['GT'] for x in record.samples]).replace('|','')).count('1')
        d_af[i] = record.INFO['AC'][0]
    l_af = [d_af[x] for x in d_af.keys()]
    return(l_af)    


#Function to extract sfs from dictionary of allele frequencies
def get_sfs(l_af, samples):
    d_sfs = {}
    s_seg = 0 #total number of truly segregating sites
    s_not_anc = 0 #required to know the d0_0 class
    #Loop through list of allele frequency categories
    for x in l_af:
        try:
            #if the category already exists, incrememnt by one
            d_sfs[x] = d_sfs[x] + 1
        except:
            #otherwise create the category
            d_sfs[x] = 1
        #add to count of segregating sites if not fixed or lost
        if int(x) > 0 and int(x) < int(samples):
            s_seg += 1
        if int(x) > 0:
            s_not_anc += 1
    return(d_sfs, s_seg, s_not_anc)


#Master function to produce dadi spectrum and output to file
def master_func(vcf_file, samples, outFile=''):
    if(sample_type=='chimeric'):
        l_af = get_chimeric_counts(vcf_file)
    else:
        l_af = get_counts(vcf_file)
    sfs_list = get_sfs(l_af, samples)
    sfs = sfs_list[0]
    sfs = dict(sorted(sfs.items()))
    s_seg = sfs_list[1]#Number of polymorphic sites, excludes sites with AF=1
    s_not_anc = sfs_list[2] #Number of polymorphic sites + number of fixed sites in the sample      
    sfs_arr = list(sfs.values())
    sfs_arr = [0] + sfs_arr
    dadi.Spectrum(sfs_arr).to_file(outFile) 
    print ("sfs output to " + outFile)


master_func(vcf_file, samples, outFile)

