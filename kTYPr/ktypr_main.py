import pandas as pd
import os, sys
import glob
import pickle
import ktypr_hmms as kh
import ktypr_utils as ku
from pathlib import Path
from joblib import Parallel, delayed
from datetime import datetime

### GLOBAL VARIABLES
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

hmms_path        = os.path.join(SCRIPT_DIR, 'data', 'hmms', '*.hmm')
kpsc_hmm_path    = os.path.join(SCRIPT_DIR, 'data', 'hmms', 'KpsC.hmm')
definition_path  = os.path.join(SCRIPT_DIR, 'data', 'ktypr_definitions_v20250512.tsv')
cutoffs_path     = os.path.join(SCRIPT_DIR, 'data', 'hmm_cutoffs_v20250704.tsv')

thrs_dict = kh.load_bitscore_thrs(cutoffs_path)
ktype_dict, subject_to_ktypes, ktype_gene_counts = kh.get_ktypes_dicts(definition_path)
kanalysis_columns, korder = kh.get_ktypr_columns_and_order(ktype_dict)

### MAIN

def set_directory_outputs(identifier, outDir, flanking=False):
    """
    Sets directories and defines paths to store hits and prediction files 
    """
    if outDir[-1]!='/':
        outDir += '/'
    outDir += f'{identifier}/'
    Path(outDir).mkdir(parents=True, exist_ok=True)
    outdict = {'outdir':outDir, 'id':identifier, 
               'hits':f'{outDir}{identifier}_hits.tsv.gz', 
               'classification':f'{outDir}{identifier}_ktypr.tsv'}
    if flanking:
        outdict['flank'] = f'{outDir}{identifier}_flanks.faa'

    return outdict

def profile_genome_from_aa_annotations(aa_file, 
                                       hmms_path=hmms_path, 
                                       thrs=thrs_dict, 
                                       conserved = ['KpsEDCSMT', 'KpsFU'], 
                                       identifier=None,
                                       outDir=None,
                                       append=None,
                                       _rfbBDAC=True):
    """
    Runs the default pipeline given a fasta aa file 
    """

    # Run HMMs
    hits = kh.retrieve_hits(aa_file, hmms_path)

    # Set directories and identifier to save the output
    if not identifier:
        identifier = aa_file.split('/')[-1].split('.')[0]
    if outDir:
        # Get identifier and set directory
        manager_dict = set_directory_outputs(identifier, outDir)
        hits.to_csv(manager_dict['hits'], sep='\t', compression='gzip')   # Save all hits

    # Filter by bitscore, group by unique subject hit and produce the dictionary with K-type info
    hits = kh.calculate_hits_and_bitscores(kh.filter_max_bitscore(kh.apply_bitscore_thresholds(hits, thrs)), 
                                           subject_to_ktypes, ktype_gene_counts, _rfbBDAC)
    
    # Get best and add it to the dictionary
    best = kh.get_best_k({k:v for k, v in hits.items() if k not in conserved})
    hits['pred'] = hits.get(best, [0, 0, 0, 0])

    # Make output
    if outDir:
        results = [identifier, best]+hits['pred']
        for k in korder:
            results += hits.get(k, [ktype_gene_counts[k],0,0,0])
        if append:
            with open(append, 'a') as fo:
                fo.write('\t'.join([str(i) for i in results])+'\n')
        else:
            with open(manager_dict['classification'], 'w') as fo:
                fo.write('\t'.join(kanalysis_columns)+'\n')
                fo.write('\t'.join([str(i) for i in results])+'\n')

    return hits


def profile_genome_from_aa_annotations_flanking(aa_file, 
                                                flank = 30000,
                                                gff_file=None,
                                                flank_hmm_path=kpsc_hmm_path, 
                                                hmms_path=hmms_path, 
                                                thrs=thrs_dict, 
                                                conserved = ['KpsEDCSMT', 'KpsFU'], 
                                                identifier=None,
                                                outDir='./',
                                                append=None,
                                                _multi_kps=False, 
                                                _rfbBDAC=False):
    """
    Runs the default pipeline given a fasta aa file 
    """

    # Check if KpsC present
    kpsc = kh.profile_KpsC(aa_file, flank_hmm_path, multi=_multi_kps)
    # Subset flanking
    if not identifier:
        identifier = aa_file.split('/')[-1].split('.')[0]

    if kpsc:
        manager_dict = set_directory_outputs(identifier, outDir, flanking=True)

        # Filter the fasta based on annotation
        ku.subset_flanking(aa_file, kpsc, gff_file=gff_file, flank=30000, out_file=manager_dict['flank'])

        # Run HMMs
        hits = kh.retrieve_hits(manager_dict['flank'], hmms_path)

        if outDir:
            # Get identifier and set directory
            hits.to_csv(manager_dict['hits'], sep='\t', compression='gzip')   # Save all hits

        # Filter by bitscore, group by unique subject hit and produce the dictionary with K-type info
        hits = kh.calculate_hits_and_bitscores(kh.filter_max_bitscore(kh.apply_bitscore_thresholds(hits, thrs)), 
                                               subject_to_ktypes, ktype_gene_counts, _rfbBDAC)
    
        # Get best and add it to the dictionary
        best = kh.get_best_k({k:v for k, v in hits.items() if k not in conserved})
    else:
        best = 'No KpsC identified'
        hits = {}

    hits['pred'] = hits.get(best, [0, 0, 0, 0])

    # Make output
    if outDir:
        results = [identifier, best]+hits['pred']
        for k in korder:
            results += hits.get(k, [ktype_gene_counts[k],0,0,0])
        if append:
            with open(append, 'a') as fo:
                fo.write('\t'.join([str(i) for i in results])+'\n')
        else:
            with open(manager_dict['classification'], 'w') as fo:
                fo.write('\t'.join(kanalysis_columns)+'\n')
                fo.write('\t'.join([str(i) for i in results])+'\n')

    return hits


# FOR MULTIPLE GENOMES

def run_on_collection(inFile,
                      outDir, 
                      collection_id, collection_folder=None,
                      flanking=False, flank=30000, multi=False, 
                      _rfbBDAC=False, 
                      parallel=True, n_jobs=10, redo=1, test=False):
    
    # Process all input files in parallel
    with open(inFile, 'r') as file:
        fils = [line.strip() for line in file]    
    
    # For testing purposes
    if test:
        fils = fils[:10]

    # Manage file to append and redo
    if collection_folder:
        # Generate the current timestamp
        timestamp = datetime.now().strftime("%Y%m%d")
        append_fil = f'{collection_folder}/{collection_id}_{timestamp}_ktypr.tsv'
    else:
        append_fil = f'{outDir}/{collection_id}_ktypr.tsv'

    if redo:
        if os.path.isfile(append_fil):
            os.remove(append_fil)
        with open(append_fil, 'w') as fo:
            fo.write('\t'.join(['genome_id']+kanalysis_columns)+'\n')
    else:
        pass
    
    # Run in parallel
    if parallel:
        if flanking:
            Parallel(n_jobs=n_jobs)(delayed(profile_genome_from_aa_annotations_flanking)(aa_file=fil, 
                                                                                         flank=flank,
                                                                                         gff_file=fil.replace('_faa/', '_gff/').replace('.faa', '.gff'), 
                                                                                         outDir=outDir, 
                                                                                         append=append_fil,
                                                                                         _multi_kps=multi, 
                                                                                         _rfbBDAC=_rfbBDAC) for fil in fils)
        else:
            # Default
            Parallel(n_jobs=n_jobs)(delayed(profile_genome_from_aa_annotations)(aa_file=fil, 
                                                                                outDir=outDir, 
                                                                                append=append_fil, 
                                                                                _rfbBDAC=_rfbBDAC) for fil in fils)

def ktypr_main(inFile, outDir, collection_id, collection_folder=None, 
               flanking=False, flank=30000, multi=False, _rfbBDAC=False, 
               parallel=True, n_jobs=10, redo=1, test=False):
    """
    Main function to run ktypr on a collection of genomes
    """
    run_on_collection(inFile, outDir, collection_id, collection_folder,
                      flanking=flanking, flank=flank, multi=multi,
                      _rfbBDAC=_rfbBDAC,
                      parallel=parallel, n_jobs=n_jobs, redo=redo, test=test)

if __name__ == "__main__":

    fils = glob.glob('/nfs/home/smiravet/KTYPS_DEV/scratch/raw/*/fetch_annotations.txt')

    torun = ['gladstone', 'fastkaptive', 'kref10', 'bioinforef', 'refseq37', 'ncbi', 'egli', 'barn', 'kref']
    torun = ['sands', 'nuin']
    torun = ['galardini', 'botnar', 'colibafi', 'coliville', 'ecoref', 'lbc', 'roar', 'septicoli', 'touchon']
    torun = ['commensal']	
    torun = ['picard']
    torun = ['hong']

    for fil in fils[::-1]:
        ide = fil.split('/raw/')[-1].split('/')[0]
        print(ide)
        if ide in torun:
            print(ide, ' running')
            if ide=='fastkaptive' or ide=='gladstone' or ide=='barn' or ide=='egli':
                # Only wg
                run_on_collection(f'/nfs/home/smiravet/KTYPS_DEV/scratch/raw/{ide}/fetch_annotations.txt', 
                                  outDir=f'/nfs/home/smiravet/KTYPS_DEV/scratch/processed/{ide}/{ide}_flank/',
                                  collection_id=f'{ide}_flank', collection_folder='/nfs/home/smiravet/KTYPS_DEV/scratch/to_report/wg_flank_v2',  n_jobs=25)
                print(ide, ' flank by default run!')
            else:
                run_on_collection(f'/nfs/home/smiravet/KTYPS_DEV/scratch/raw/{ide}/fetch_annotations.txt', 
                                  outDir=f'/nfs/home/smiravet/KTYPS_DEV/scratch/processed/{ide}/{ide}_flank/',
                                  multi=True, flanking=True, 
                                  collection_id=f'{ide}_flank', collection_folder='/nfs/home/smiravet/KTYPS_DEV/scratch/to_report/wg_flank_v2',  n_jobs=25)
                print(ide, ' flanking run!')
                run_on_collection(f'/nfs/home/smiravet/KTYPS_DEV/scratch/raw/{ide}/fetch_annotations.txt', 
                                  outDir=f'/nfs/home/smiravet/KTYPS_DEV/scratch/processed/{ide}/{ide}_wg/',
                                  collection_id=f'{ide}_wg', collection_folder='/nfs/home/smiravet/KTYPS_DEV/scratch/to_report/wg_flank_v2',  n_jobs=25)
                print(ide, ' wg run!')
