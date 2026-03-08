import pandas as pd
import os, sys
import glob
import pickle
from . import ktypr_hmms as kh
from . import ktypr_utils as ku
from pathlib import Path
from joblib import Parallel, delayed
from datetime import datetime

### GLOBAL VARIABLES
SCRIPT_DIR        = os.path.dirname(os.path.abspath(__file__))
hmms_path         = os.path.join(SCRIPT_DIR, 'data', 'hmms', '*.hmm')
kpsc_hmm_path     = os.path.join(SCRIPT_DIR, 'data', 'hmms', 'KpsC.hmm')
definition_path   = os.path.join(SCRIPT_DIR, 'data', 'ktypr_definitions_v20260308.tsv')
cutoffs_path      = os.path.join(SCRIPT_DIR, 'data', 'hmm_cutoffs_v20260308.tsv')
thrs_dict         = kh.load_bitscore_thrs(cutoffs_path)
max_bitscore_dict = kh.load_hmm_bitscore_max()    # Path to the file defined as default in the function

# Load definitions and report info
ktype_dict, subject_to_ktypes, ktype_gene_counts = kh.get_ktypes_dicts(definition_path)
kanalysis_columns, korder                        = kh.get_ktypr_columns_and_order(ktype_dict)

# Profiling functions

def _profile_genome_core(faa_path,
                         result_dict,
                         hmms_path,
                         thrs,
                         conserved,
                         append,
                         _rfbBDAC,
                         classification_label='',
                         skip_hits=False, 
                         ignore_cutoffs=False):
    """
    Core function to run HMM search, process hits, and write classification.
    """
    if skip_hits:
        best = classification_label or 'No KpsC identified'
        hits = {}
        max_hits = None
    else:
        # Run HMMs
        hits_df = kh.retrieve_hits(faa_path, hmms_path)
        hits_df.to_csv(result_dict['hits'], sep='\t', compression='gzip')  # Save all hits

        # Filter by bitscore, group by unique subject hit and produce the dictionary with K-type info
        
        if ignore_cutoffs:
            print("WARNING: Ignoring bitscore cutoffs for HMM prediction. This may lead to more false positives. Use with caution.")
            thrs = {k: 0 for k in thrs}  # Set all thresholds to 0 to effectively ignore them

        filtered_hits = kh.apply_bitscore_thresholds(hits_df, thrs)
        max_hits = kh.filter_max_bitscore(filtered_hits)
        
        # Save filt_hits
        max_hits.to_csv(result_dict['filt_hits'], sep='\t', compression='gzip')  # Save filtered hits

        # Count hits per ktype
        hits = kh.calculate_hits_and_bitscores(df=max_hits, subject_to_ktypes=subject_to_ktypes, 
                                               ktype_gene_counts=ktype_gene_counts, max_bitscore_dict=max_bitscore_dict,
                                               _rfbBDAC=False)

        # Get best prediction
        best = kh.get_best_k({k: v for k, v in hits.items() if k not in conserved})

    # Store predicted best result
    hits['pred'] = hits.get(best, [0, 0, 0, 0, 0])

    # Store output
    results = [result_dict['ide'], best] + hits['pred']
    for k in korder:
        results += hits.get(k, [ktype_gene_counts[k], 0, 0, 0, 0])

    with open(result_dict['classification'], 'w') as fo:
        fo.write('\t'.join(kanalysis_columns) + '\n')
        fo.write('\t'.join(str(i) for i in results) + '\n')

    if append:
        with open(append, 'a') as fo:
            fo.write('\t'.join(str(i) for i in results) + '\n')

    # Update
    result_dict['max_hits'] = max_hits
    result_dict['done'] = 1
    result_dict['ktype'] = best

    return hits

def profile_genome_from_aa_annotations(result_dict, 
                                       hmms_path=hmms_path, 
                                       thrs=thrs_dict, 
                                       conserved=['KpsEDCSMT', 'KpsFU'], 
                                       append=None,
                                       _rfbBDAC=True, 
                                       ignore_cutoffs=False):
    """
    Runs the default pipeline given a fasta aa file.
    """
    return _profile_genome_core(result_dict['faa_path'],
                                result_dict,
                                hmms_path,
                                thrs,
                                conserved,
                                append,
                                _rfbBDAC,
                                ignore_cutoffs=ignore_cutoffs)


def profile_genome_from_aa_annotations_flanking(result_dict, 
                                                flank=30000,
                                                flank_hmm_path=kpsc_hmm_path, 
                                                hmms_path=hmms_path, 
                                                thrs=thrs_dict, 
                                                conserved=['KpsEDCSMT', 'KpsFU'], 
                                                append=None,
                                                _multi_kps=False, 
                                                _rfbBDAC=False,
                                                ignore_cutoffs=False):
    """
    Runs the default pipeline only for annotations flanking the KpsC gene.
    """
    kpsc = kh.profile_KpsC(result_dict['faa_path'], flank_hmm_path, multi=_multi_kps)
    if kpsc:
        ku.subset_flanking(result_dict['faa_path'], kpsc,
                           gff_file=result_dict['gff_path'],
                           flank=flank,
                           out_file=result_dict['faa_flank_path'])

        return _profile_genome_core(result_dict['faa_flank_path'],
                                    result_dict,
                                    hmms_path,
                                    thrs,
                                    conserved,
                                    append,
                                    _rfbBDAC,
                                    ignore_cutoffs=ignore_cutoffs)
    else:
        return _profile_genome_core(result_dict['faa_path'],
                                    result_dict,
                                    hmms_path,
                                    thrs,
                                    conserved,
                                    append,
                                    _rfbBDAC,
                                    classification_label='No KpsC identified',
                                    skip_hits=True, 
                                    ignore_cutoffs=ignore_cutoffs
                                    )


def annotate_and_profile(genome_path, outDir, prefix='', extract_annotations=True, 
                         reannotate=False, meta=False, flanking=False, flank=30000,
                         multi=False, _rfbBDAC=False, append_fil=None, clinker=True, 
                         verbose=True, ignore_cutoffs=False):
    """
    Annotates a genome or annotation file and profiles it for K-type prediction.
    """    
    result = ku.prepare_single_input(
        genome_path=genome_path, 
        outDir=outDir, 
        prefix=prefix, 
        extract_annotations=extract_annotations, 
        reannotate=reannotate, 
        meta=meta,
        verbose=verbose
    )

    if flanking:
        if verbose:
            print('Profiling genome with flanking KpsC annotations...')
        profile_genome_from_aa_annotations_flanking(
            result_dict=result,
            flank=flank,
            append=append_fil,
            _multi_kps=multi,
            _rfbBDAC=_rfbBDAC, 
            ignore_cutoffs=ignore_cutoffs
        )
    else:
        result['faa_flank_path'] = None
        if verbose:
            print('Profiling genome with whole-genome annotations...')
        profile_genome_from_aa_annotations(
            result_dict=result,
            append=append_fil,
            _rfbBDAC=_rfbBDAC,
            ignore_cutoffs=ignore_cutoffs
        )

    # Report genbank and clinker
    if result['gff_path'] is not None:
        if verbose:
            print('Creating GenBank file...')
        ku.create_genbank_from_inputs(result_dict=result,
            id_attribute="ID",
            from_annotations=False, 
            verbose=verbose)
                
    return result

# MAIN FUNCTION

def ktypr(inFile, outDir, prefix='',
          flanking=False, flank=30000, 
          reannotate=False, multi=False, _rfbBDAC=False, 
          parallel=True, n_jobs=10, 
          meta=False, clinker=True, 
          verbose=False, keep_output=True, ignore_cutoffs=False):
    """
    Main function to run ktypr on one or multiple genomes.

    Parameters
    ----------
    inFile : str or Path or list
        Input FASTA or GENBANK file(s), list with paths, or directory containing genome sequences to process.
    outDir : str or Path
        Output directory where results and intermediate files will be saved.
    prefix : str, optional
        Prefix for output file names. If multiple genomes are processed, this prefix is 
        prepended to resulting files.
    flanking : bool, default=False
        Whether to include and analyze flanking genes (+-flank) around detected loci.
    flank : int, default=30000
        Number of base pairs to include upstream and downstream when flanking=True.
    reannotate : bool, default=False
        If True, forces re-annotation of genomes even if previous annotations are found such as in a full genbank file.
    multi : bool, default=False
        If True, allows for finding more than one cluster in flanking mode, for example if multiple kpsC hits are found. 
    _rfbBDAC : bool, default=False
        Experimental option to account for the scoring genes rfbB, D, A, and C.
    parallel : bool, default=True
        Whether to run genome profiling in parallel using joblib when multiple genomes are provided.
    n_jobs : int, default=10
        Number of parallel workers to use if parallel=True.
    meta : bool, default=False
        Whether the input genomes are metagenome-assembled (MAGs) or incomplete drafts, important in the annotation step as
        it will set meta=True for the pyrodigal gene calling. 
    clinker: bool, default=True
        Runs a clinker analysis with the produced genbank against the predicted k-type.
    verbose : bool, default=False
        If True, prints progress and intermediate information to stdout.
    keep_output : bool, default=True
        If True, keeps intermediate output files. If False, removes them after processing.
    ignore_cutoffs : bool, default=False
        If True, ignores cutoffs for HMM prediction. Use with caution, as this may lead to more false positives.

    Output
    -----
    One folder per genome including annotations, clinker

    Notes
    -----
    This function performs annotation and K-type profiling on a collection of genome 
    sequences. When multiple input files are provided, results are appended to a single
    TSV file. Parallel processing is supported via joblib.
    """

    # Resolve input file(s) and check they are valid
    input_paths = ku.resolve_paths(inFile, verbose=verbose)
    if verbose:
        print(f"Resolved input paths: {input_paths}")

    if len(input_paths) > 1:
        os.makedirs(outDir, exist_ok=True)
        if prefix:
            append_fil = os.path.join(outDir, f'{prefix}results_ktypr.tsv')
        else:
            append_fil = os.path.join(outDir, 'results_ktypr.tsv')
            prefix = ''  # default to empty string if not given

        # Write header to the file
        with open(append_fil, 'w') as fo:
            fo.write('\t'.join(kanalysis_columns) + '\n')
    else:
        append_fil = None

    # Annotate and profile each genome
    if parallel and n_jobs!=1:
        results = Parallel(n_jobs=n_jobs, backend="multiprocessing")(
                      delayed(annotate_and_profile)(
                          genome_path=genome_path,
                          outDir=outDir,
                          prefix=prefix,
                          extract_annotations=True,
                          reannotate=reannotate,
                          meta=meta,
                          flanking=flanking,
                          flank=flank,
                          multi=multi,
                          _rfbBDAC=_rfbBDAC,
                          append_fil=append_fil,
                          verbose=verbose, 
                          ignore_cutoffs=ignore_cutoffs
                              ) for genome_path in input_paths
                      )
    else:
        # This is only when no parallel or n_jobs=1
        # Not really needed, but should be faster than n_jobs=1 as it skips joblib overhead
        results = []
        for genome_path in input_paths:
            res = annotate_and_profile(
                          genome_path=genome_path,
                          outDir=outDir,
                          prefix=prefix,
                          extract_annotations=True,
                          reannotate=reannotate,
                          meta=meta,
                          flanking=flanking,
                          flank=flank,
                          multi=multi,
                          _rfbBDAC=_rfbBDAC,
                          append_fil=append_fil,
                          verbose=verbose,
                          ignore_cutoffs=ignore_cutoffs
                     )   
        results.append(res)

    # Clinker collectively
    if clinker and flanking:
        ku.get_clinker(results, verbose=verbose)
        
    # Cleanup intermediate files if needed
    if keep_output:
        return results
    else:
        ku.cleanup_intermediate_files(results, verbose=verbose)
        results.clear()
