# Repository for Soni et al. 2025: Investigating the effects of chimerism on the inference of selection: quantifying genomic targets of purifying, positive, and balancing selection in common marmosets (Callithrix jacchus)
Vivak Soni, Cyril J. Versoza, Susanne P. Pfeifer, Jeffrey D. Jensen


empirical_scans: Contains input files for running selective sweep scans with SweepFinder2 and balancing selection scans with B0MAF, as well as scan results. Also contains mean gene divergence file.


power_analysis: Contains results files from selective sweep and balancing selection scans on simulated data of an equilibrium population (WF and chimerism), and the marmoset demographic model (chimerism only). 


DFE_inference.tar.gz: Contains summary statistics from 100 simulation replicates of best-fitting DFE under marmoset demographic model, as well as fixation files to calculate divergence. Also contains file of empirical exonic divergence.


scripts: Contains SLiM and python scripts to perform all analyses:
    
    marm_bestFit.slim: SLiM script to simulate the marmoset demographic model inferred by Soni et al. 2025 (MBE).
    
    marm_bestFit_bs.slim: SLiM script to simulate marmoset demographic model under balancing selection.
    
    marm_bestFit_DFE.slim: SLiM script to simulate marmoset demographic model with divergence phase for DFE inference.
    
    marm_bestFit_sweeps.slim: SLiM script to simulate marmoset demographic model under a selective sweep.
    
    WF_bs.slim: SLiM script to simulate WF population under balancing selection.
    
    WF_sweeps.slim: SLiM script to simulate WF population under a selective sweep.
    
    get_samples.py: Python script to get samples of chimeras (see https://github.com/vivaksoni/marmoset_chimerism_demography for further information)
    
    get_input_for_BM.py: Python script to generate input file for B0MAF balancing selection scans.
    
    get_input_for_SF2.py: Python script to generate input file for SweepFinder2 selective sweep scans.
    
    get_summary_stats.py: Python script to generate summary statistics across sliding windows from vcf input file.
    
    marmoset_selection_plots.ipynb: Jupyter notebook for generating main figures from manuscript.
    
    slim4_chimerism_DFE.sh and .sub: Example submission files for submitting jobs on Open Science Grid, showing pipeline for running simulations and obtaining chimeric samples.
