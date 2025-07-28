import os
import glob
import gzip
import pyhmmer
import collections
import pandas as pd
import ktypr_hmms as kh
from itertools import chain
from pyhmmer.easel import SequenceFile, Alphabet

### GLOBAL VARIABLES
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

protein_alphabet = Alphabet.amino()
definition_path  = os.path.join(SCRIPT_DIR, 'data', 'ktypr_definitions_v20250512.tsv')
cutoffs_path     = os.path.join(SCRIPT_DIR, 'data', 'hmm_cutoffs_v20250704.tsv')
kpsc_hmm_path    = os.path.join(SCRIPT_DIR, 'data', 'hmms', 'KpsC.hmm')

#definition_path = '/nfs/home/smiravet/KTYPS_DEV/data/kTYP_definitions_v241202.tsv'
#cutoffs_path = '/nfs/home/smiravet/KTYPS_DEV/data/hmm_cutoffs_v241201.tsv'

### I/O

def load_k_assignments(inFile=definition_path):
    """
    Given a text file <inFile> in tsv format including which HMMs relate to each KTYPE such as:
    K1\thmm1_1\thmm1_2\t...\thmm1_N\n
    K2\thmm2_1\thmm2_2\t...\thmm2_N\n

    Returns a dictionary {K1:[hmm1_1, hmm1_2], etc.}

    Notice the same hmm can contribute to different Ktypes
    """
    assignments = {}
    with open(inFile, 'r') as fi:
        for line in fi:
            line = line.strip().split('\t')
            assignments[line[0]] = line[1:]
    return assignments


def get_ktypes_dicts(inFile=definition_path):

    ktype_dict = load_k_assignments(inFile)
    # Reverse mapping: create a dictionary where each subject points to all relevant ktypes
    subject_to_ktypes = collections.defaultdict(list)
    ktype_gene_counts = {}
    for ktype, subjects in ktype_dict.items():
        ktype_gene_counts[ktype] = len(subjects)
        for subject in subjects:
            subject_to_ktypes[subject].append(ktype)
    
    return ktype_dict, subject_to_ktypes, ktype_gene_counts


def load_bitscore_thrs(inFile=cutoffs_path):
    """
    Given a text file <inFile> in tsv format including which HMMs and given bitscore thresholds:
    K1\tthr\n
    K2\tthr\n
    """
    thrs = {}
    header = True
    with open(inFile, 'r') as fi:
        for line in fi:
            if header:
                header = False
            else:
                k, v = line.strip().split()
                thrs[k] = int(v)
    return thrs


### HMMs

def retrieve_hits(seqs_path, hmms_path, 
                  fields=["query", "subject", "bitscore", "evalue"],  max_aa_size=2000):
    """
    Given a sequence path <seqs_path> in fasta format and a directory with hmms <hmms_path>,
    returns a dataframe with the hits containing <fields>
    """
    # Load cluster proteins
    with pyhmmer.easel.SequenceFile(seqs_path, digital=True, alphabet=protein_alphabet) as seqs_file:
        proteins = [seq for seq in list(seqs_file) if len(seq) <= max_aa_size] 

    # Load HMMs
    hmms = []
    for fil in glob.glob(hmms_path):
        with pyhmmer.plan7.HMMFile(fil) as hmm_file:
            hmms.append(hmm_file.read())

    # Run HMMs
    Result = collections.namedtuple("Result", fields)

    results = []
    for hits in pyhmmer.hmmsearch(hmms, proteins, E=1):
        cog = hits.query_name.decode()
        for hit in hits:
            if hit.included:
                results.append(Result(hit.name.decode(), cog, hit.score, hit.evalue))

    # Results --> df
    hits_df = {}
    c = 0
    for i in results:
        hits_df[c] = list(i)
        c += 1
    hits_df = pd.DataFrame.from_dict(hits_df, orient='index', columns=fields)

    return hits_df


def apply_bitscore_thresholds(df, thrs, _rfbBDAC=True):
    """ Filters out results not passing the cut-offs. """
    filter_col = []
    for _, row in df.iterrows():
        if row['bitscore']>thrs.get(row['subject'], 0):
            filter_col.append(1)
        else:
            filter_col.append(0)
    df['bitscore_thr'] = filter_col
    filt_df = df[df['bitscore_thr']==1]
    return filt_df


def filter_max_bitscore(df, group_by='subject', by='bitscore'):
    """ Filter the DataFrame to keep only the rows with maximum bitscore per query """
    max_bitscore_idx = df.groupby(group_by)[by].idxmax()   # Find max
    filtered_df = df.loc[max_bitscore_idx].reset_index(drop=True) 
    return filtered_df


def calculate_hits_and_bitscores(df, subject_to_ktypes, ktype_gene_counts, _rfbBDAC=False):
    """
    Example:
    subject_to_ktypes = {'K1a': ['K1', 'K2'], 'K1b': ['K1'], 'K2a': ['K2'], 'RfbB':['K1']}
    ktype_gene_counts = {'K1':3, 'K2':2}
    df = pd.DataFrame({
        'subject': ['K1a', 'K1a', 'K1b', 'K2a', 'RfbB'],
        'bitscore': [200, 50, 100, 50, 1000] 
    })
    calculate_hits_and_bitscores(df, subject_to_ktypes, ktype_gene_counts, False) --> {'K1': [3, 4, 350, 0], 'K2': [2, 3, 300, 0]}

    If the df is filtered for unique subjects, the second value in the return dictionary will always be <= than number of genes in
    the ktype. 
    """

    # Initialize dictionaries for counting hits and accumulating bitscores
    hit_counts = collections.defaultdict(int)
    accumulated_bitscores = collections.defaultdict(int)
    
    # Iterate over each row in the DataFrame
    rs = {}
    for _, row in df.iterrows():
        subject, bitscore = row['subject'], row['bitscore']
        for ktype in subject_to_ktypes[subject]:
            if ktype not in rs:
                rs[ktype] = [ktype_gene_counts[ktype], 0, 0, 0]   # Start the results counter
            # Update hit count and bitscore for the corresponding ktype
            rs[ktype][1] += 1   # Add gene to count of gene, rfb genes should count for completeness but not for bitscore 
            if _rfbBDAC:
                rs[ktype][2] += bitscore  # Add up bitscore independently on the hit
            else:
                if subject not in {'RfbB', 'RfbD', 'RfbA', 'RfbC'}:
                    rs[ktype][2] += bitscore   # Add up only if no Rfb gene (still count them for the ktype number of genes, but the bitscore doesn't)
    # Add complete
    for k, v in rs.items():
        if v[1] >= v[0]:
            v[-1] = 1
    return rs
    

def get_best_k(ktype_analysis, not_found_default='No K-type assigned'):
    """
    Given a dictionary of ktype_analysis such as {'K1': [3, 3, 350, 0], 'K2': [2, 1, 300, 0]}
    add to the dictionary a 'best' entry with the maximum bitscore k-type that is complete (if no complete just the max bitscore)
    if add_to_dict false, just return the Ktype

    If empty input, reports 'No K-type assigned' 
    Be sure the conserved regions are not included here!
    """
    
    if not ktype_analysis:
        return not_found_default
    # Use max with a custom key to find the best ktype
    best_ktype, _ = max(
        ktype_analysis.items(),
        key=lambda item: (item[1][3], item[1][2]),  # Prioritize by completeness (index 3), then bitscore (index 2)
        default=(not_found_default, None)  # Handle empty input gracefully
    )
    return best_ktype


def get_multiple_best_k(ktype_analysis):
    """
    Determines the best k-types based on completeness and, if necessary,
    the number of genes present in the genome. Allows multiple k-types 
    if they tie on the criteria.
    
    Parameters:
        ktype_analysis (dict): A dictionary with k-type keys and lists as values.
            Each list contains:
            [nr_of_genes_in_ktype, nr_of_genes_from_ktype_in_genome, accumulated_bitscore, is_complete].
    
    Returns:
        dict: A dictionary of the best k-types and their corresponding data.
    
    Be sure the conserved regions are not included here!
    """
    if not ktype_analysis:
        return {'No K-type assigned': None}
    
    # Group k-types by completeness
    complete_ktypes = {
        ktype: data for ktype, data in ktype_analysis.items() if data[3] > 0
    }
    
    # If there are complete k-types, filter for the maximum completeness
    if complete_ktypes:
        max_completeness = max(data[3] for data in complete_ktypes.values())
        best = {
            ktype: data for ktype, data in complete_ktypes.items() if data[3] == max_completeness
        }
    else:
        # Otherwise, use the number of genes present as the criterion
        max_genes_present = max(
            data[1] for data in ktype_analysis.values()
        )
        best = {
            ktype: data for ktype, data in ktype_analysis.items() if data[1] == max_genes_present
        }
    
    return best


### Specific functions for flanking profiling

def profile_KpsC(seqs_path, hmm_path=kpsc_hmm_path, multi=False):
    """ Profiles for KpsC used as marker gene for the cluster presence """
    df = retrieve_hits(seqs_path, hmm_path)
    if df.empty:
        return False
    else:
        if multi:
            return list(df['query'])
        else:
            return list([df.loc[df['bitscore'].idxmax(), 'query']])

### Analysis and report

def get_ktypr_columns_and_order(assignments):
    """
    Get columns for the classification table
    """
    # Define conserved types and dynamic types (ktype) based on assignments
    avail_cons = ['KpsEDCSMT', 'KpsFU']
    avail_ktyp = sorted(set(assignments) - set(avail_cons))
    
    # Define common suffixes for each column
    suffixes = ['_nr_genes', '_genes_in_genome', '_acc_bitscore', '_is_complete']
    
    # Construct columns list using comprehensions
    columns = ['predicted', 'pred_nr_genes', 'pred_genes_in_genome', 'pred_acc_bitscore', 'pred_is_complete']
    columns += [f"{cons}{suffix}" for cons in avail_cons for suffix in suffixes]
    columns += [f"{ktyp}{suffix}" for ktyp in avail_ktyp for suffix in suffixes]
    
    return columns, avail_cons+avail_ktyp