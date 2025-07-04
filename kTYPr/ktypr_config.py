# code/config.py
from pathlib import Path
import ktyps_hmms as kh

ROOT_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = ROOT_DIR / 'data'

HMM_CUTOFFS_PATH = DATA_DIR / 'hmm_cutoffs_v241201.tsv'
KTYP_DEFS_PATH = DATA_DIR / 'kTYP_definitions_v241202.tsv'
HMMS_PATH = DATA_DIR / 'hmms' / '*.hmm'
KPSC_HMM_PATH = DATA_DIR / 'hmms' / 'KpsC.hmm'

# Load once and import as needed
thrs_dict = kh.load_bitscore_thrs(HMM_CUTOFFS_PATH)
ktype_dict, subject_to_ktypes, ktype_gene_counts = kh.get_ktypes_dicts(KTYP_DEFS_PATH)
kanalysis_columns, korder = kh.get_ktyps_columns_and_order(ktype_dict)