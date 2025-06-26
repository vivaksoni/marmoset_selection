#!/bin/bash

# have job exit if any command returns with non-zero exit status (aka failure)
set -e

# replace env-name on the right hand side of this line with the name of your conda environment
ENVNAME=env5

# if you need the environment directory to be named something other than the environment name, change this line
ENVDIR=$ENVNAME

# these lines handle setting up the environment; you shouldn't have to modify them
export PATH
mkdir $ENVDIR
tar -xzf $ENVNAME.tar.gz -C $ENVDIR
. $ENVDIR/bin/activate

# modify environment variables
export PATH=$_CONDOR_SCRATCH_DIR/build:$PATH

tar -zxf scripts.tar.gz

slim -d "d_paramID='$1'" -d "d_repID='$2'" -d "d_folder='.'" -d d_f0=$3 -d d_f1=$4 -d d_f2=$5 -d d_f3=$6 -d d_mu=$7 -d d_r=$8 -d d_L=$9 scripts/marm_bestFit_DFE.slim

python3 scripts/get_samples.py -pedigrees pedigreeIDs_$1_rep$2.txt -outFile $1_rep$2 -samples 15

vcftools --vcf $1_rep$2.vcf --keep $1_rep$2_chimeric.txt --out $1_rep$2_chimeric --recode

rm $1_rep$2.vcf
mv $1_rep$2_chimeric.recode.vcf $1_rep$2_chimeric.vcf

python3 scripts/get_summary_stats_noWins.py -inFile $1_rep$2_chimeric.vcf -outFile $1_rep$2_chimeric.stats -regionLen 3209

mkdir DFE_${10}
cp *.txt DFE_${10}
cp *.vcf DFE_${10}
cp *.stats DFE_${10}
cp *.fixed DFE_${10}
tar -czf DFE_${10}.tar.gz DFE_${10}
