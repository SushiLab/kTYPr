import os
import glob
import gzip
import pyrodigal
import pandas as pd
import ktypr_hmms as kh
from joblib import Parallel, delayed
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import FeatureLocation
from pathlib import Path
import zipfile

### GLOBAL VARIABLES
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

### I/O

def set_pyrodigal_model(records, meta=False):
    """ Sets up the Pyrodigal gene finder model based on the provided records."""
    if meta:
        orf_finder = pyrodigal.GeneFinder(meta=True)
    else:
        full_genome_seq = b''.join([bytes(record.seq) for record in records])  # Concatenate all sequences into one byte string
        orf_finder = pyrodigal.GeneFinder()
        orf_finder.train(full_genome_seq)
    return orf_finder

def get_faa_and_gff(fil, outDir_faa, outDir_gff):
    """ Runs prodigal on a genome file. It accepts both compressed (as gz) and uncompressed fasta files """

    ide = fil.split('/')[-1].split('.')[0]
    faa_basename = f'{outDir_faa}/{ide}.faa'
    gff_basename = f'{outDir_gff}/{ide}.gff'
    ffo = open(faa_basename, "w")
    gfo = open(gff_basename, "w")

    if fil.endswith('.gz'):
        with gzip.open(fil, "rt") as fil:
            records = SeqIO.parse(fil, "fasta")
            orf_finder = set_pyrodigal_model(records, meta=False)
            for i, record in enumerate(records):
                genes = orf_finder.find_genes(bytes(record.seq))
                genes.write_gff(gfo, sequence_id=f'{ide}__{record.id}', header=(i == 0))     # Do not include header 
                genes.write_translations(ffo, sequence_id=f'{ide}__{record.id}')
    else:
        records = SeqIO.parse(fil, "fasta")
        orf_finder = set_pyrodigal_model(records, meta=False)
        for i, record in enumerate(records):
            genes = orf_finder.find_genes(bytes(record.seq))
            genes.write_gff(gfo, sequence_id=f'{ide}__{record.id}', header=(i == 0))     # Do not include header 
            genes.write_translations(ffo, sequence_id=f'{ide}__{record.id}')
    
    ffo.close()
    gfo.close()


def get_ann_dict_from_fasta(fasta_file):
    return {
        parts[0]: [int(i) for i in parts[1:4]]
        for record in SeqIO.parse(fasta_file, "fasta")
        if (parts := record.description.split(' # '))
    }


def get_ann_dict_from_gff(gff_file):
    df = pd.read_csv(gff_file, comment='#', sep='\t', header=None, usecols=[0, 3, 4, 6, 8])
    df[0] = df[8].str.extract(r'ID=([^;]+)')[0]
    return df.set_index(df.columns[0]).apply(lambda row: row[0:3].tolist(), axis=1).to_dict()


def subset_fasta(fasta_file, id_list, out_file=None):
    """
    Subset sequences from a FASTA file based on a given list of IDs.

    Parameters:
    - fasta_file (str): Path to the input FASTA file.
    - id_list (list): List of sequence IDs to retain.
    - output_file
    """

    filtered_seqs = (seq for seq in SeqIO.parse(fasta_file, "fasta") if seq.id in id_list)
    if out_file:
        SeqIO.write(filtered_seqs, out_file, "fasta")
    return filtered_seqs        


def subset_flanking(fasta_file, gene_id, gff_file=None, flank=30000, out_file=None):
    
    if gff_file and gff_file.endswith('gff'):
        ides = get_ann_dict_from_gff(gff_file)
    else:
        ides = get_ann_dict_from_fasta(fasta_file)
    
    flank_ides = []
    for ide in gene_id:
        st, en, _ = ides.get(ide)
        flank_start, flank_end = st - flank, en + flank
        flank_ides += [k for k, v in ides.items() if (flank_start <= v[0] <= flank_end) or (flank_start <= v[1] <= flank_end)]
    
    subset_fasta(fasta_file, set(flank_ides), out_file)

### 

def split_multigenbank_to_faa_and_gff(input_genbank,  outDir_faa, outDir_gff):
    """
    Splits a multi-GenBank file into individual genome files and extracts:
    - A `.faa` file with the CDS (coding sequences) in protein format.
    - A `.gff` file with annotations for the genes.
    
    Args:
        input_genbank (str): Path to the multi-GenBank file.
        output_dir (str): Path to the directory to store output files.
    """
    # Create output directories for faa and gff
    outDir_faa.mkdir(parents=True, exist_ok=True)
    outDir_gff.mkdir(parents=True, exist_ok=True)

    for record in SeqIO.parse(input_genbank, "genbank"):
        genome_id = record.name
        faa_file = outDir_faa / f"{genome_id}.faa"
        gff_file = outDir_gff / f"{genome_id}.gff"

        with open(faa_file, "w") as faa_out, open(gff_file, "w") as gff_out:
            # Write GFF header
            gff_out.write("##gff-version 3\n")
            gff_out.write(f"##sequence-region {genome_id} 1 {len(record)}\n")
            
            for feature in record.features:
                if feature.type == "CDS":
                    # Extract nucleotide sequence of CDS
                    cds_seq = feature.extract(record.seq)
                    
                    # Translate CDS to protein
                    protein_seq = cds_seq.translate(to_stop=True)
                    
                    # Create a protein record for the .faa file
                    protein_id = feature.qualifiers.get("locus_tag", ["unknown"])[0]
                    gene_id = feature.qualifiers.get("gene", ["unknown"])[0]

                    faa_record = SeqRecord(protein_seq, id=gene_id+'_'+protein_id, description="")
                    SeqIO.write(faa_record, faa_out, "fasta")

                    # Write GFF line
                    start = int(feature.location.start) + 1  # 1-based indexing
                    end = int(feature.location.end)
                    strand = "+" if feature.location.strand == 1 else "-"
                    attributes = f"ID={protein_id};"
                    if "product" in feature.qualifiers:
                        attributes += f"product={feature.qualifiers['product'][0]};"
                    gff_line = f"{genome_id}\tGenBank\tCDS\t{start}\t{end}\t.\t{strand}\t0\t{attributes}\n"
                    gff_out.write(gff_line)

# I/O

def process_genome_file(genome_path, outDir_faa, outDir_gff):
    genome_path = Path(genome_path)
    if genome_path.suffix in {".gb", ".gbk", ".gbff"}:
        # GenBank format
        split_multigenbank_to_faa_and_gff(str(genome_path), outDir_faa, outDir_gff)
    elif genome_path.suffix in {".fna", ".fa", ".fasta"}:
        # FASTA format
        get_faa_and_gff(str(genome_path), outDir_faa=outDir_faa, outDir_gff=outDir_gff)
    else:
        print(f"Unknown file format for {genome_path}, skipping.")

def prepare_input(inFile, outDir, prefix='', extract_annotations=True, n_jobs=30):
    """
    Prepares a list of annotation files from the given input:
    - If input is annotations (.faa), return directly.
    - If input is genome files (FASTA/GenBank), extract annotations (optionally using a user-defined prefix).
    """
    def resolve_paths(inFile):
        if inFile.endswith(".txt"):
            with open(inFile) as f:
                return [line.strip() for line in f if line.strip()]
        elif os.path.isdir(inFile):
            extensions = ("*.fasta", "*.fa", "*.fna", "*.gb", "*.gbk", "*.faa")
            files = []
            for ext in extensions:
                files.extend(glob.glob(os.path.join(inFile, ext)))
            return sorted(files)
        elif os.path.isfile(inFile):
            return [inFile]
        else:
            raise ValueError(f"Invalid input: {inFile} is neither a file nor a directory.")

    input_paths = resolve_paths(inFile)

    # If they are all .faa files, we assume annotations are already provided
    if all(path.endswith(".faa") for path in input_paths):
        return input_paths

    elif extract_annotations:
        basedir = os.path.dirname(os.path.abspath(inFile))
        ide = os.path.basename(basedir)

        fetch_txt_path = os.path.join(basedir, "fetch.txt")
        with open(fetch_txt_path, "w") as f:
            for path in input_paths:
                f.write(path + "\n")

        faa_o = os.path.join(outDir, f"{prefix}{ide}_faa")
        gff_o = os.path.join(outDir, f"{prefix}{ide}_gff")
        os.makedirs(faa_o, exist_ok=True)
        os.makedirs(gff_o, exist_ok=True)

        # Run gene calling
        Parallel(n_jobs=n_jobs)(
            delayed(process_genome_file)(genome_path, outDir_faa=faa_o, outDir_gff=gff_o)
            for genome_path in input_paths
        )

        print(f"{ide} done extracting annotations!")

        annotations_list = os.path.join(outDir, "fetch_annotations.txt")
        cmd = f'{SCRIPT_DIR}/fetch_paths.sh {faa_o} ".*\\.faa$" {annotations_list}'
        os.system(cmd)

        with open(annotations_list) as f:
            return [line.strip() for line in f if line.strip()]

    else:
        raise ValueError("Input does not appear to be annotations, and annotation extraction is disabled.")

# CLINKER 

def get_clinker(outDir, ks_to_extract):
    """
    Extract specific K-antigen clusters from a zipped archive of genbanks that can be then used with Clinker. 
    ks_to_extract is a set of substrings to search for in the filenames (eg. {"K1_", "K2_"}). 
    """

    archive_path = f"{SCRIPT_DIR}/data/clinker_clusters.zip"
    with zipfile.ZipFile(archive_path, 'r') as zipf:
        for file in zipf.namelist():
            if any(ks in file for ks in ks_to_extract):
                zipf.extract(file, path=outDir)
                print(f"Extracted: {file}")



if __name__ == "__main__":

    # FOR RUNNING PRODIGAL    
    directories = glob.glob('/nfs/home/smiravet/KTYPS_DEV/scratch/raw/fastkaptive/fetch.txt')
    directories = ['/nfs/home/smiravet/KTYPS_DEV/scratch/raw/nuin/fetch.txt', '/nfs/home/smiravet/KTYPS_DEV/scratch/raw/sands/fetch.txt']
    directories = [f'/nfs/home/smiravet/KTYPS_DEV/scratch/raw/{dataset}/fetch.txt' for dataset in ['galardini', 'botnar', 'colibafi', 'coliville', 'ecoref', 'lbc', 'roar', 'septicoli', 'touchon']]
    directories = ['/nfs/home/smiravet/KTYPS_DEV/scratch/raw/commensal/fetch.txt']
    directories = ['/nfs/home/smiravet/KTYPS_DEV/scratch/raw/picard/fetch.txt']
    directories = ['/nfs/home/smiravet/KTYPS_DEV/scratch/raw/hong/fetch.txt']

    wd = {}
    for i in directories:
        ide = i.split('/')[-2]
        basedir = i.replace('/fetch.txt', '')
        with open(i, 'r') as fi:
            faa_o, gff_o = f'{basedir}/{ide}_faa/', f'{basedir}/{ide}_gff/'
            os.makedirs(faa_o, exist_ok=True)  # succeeds even if directory exists.
            os.makedirs(gff_o, exist_ok=True)  # succeeds even if directory exists.
            wd[ide] = [[line.rstrip('\n') for line in fi], faa_o, gff_o] 
    
    # Extract annotation for all collections   
    for ide, info in wd.items():
        print(ide, len(info[0]))
        Parallel(n_jobs=30)(delayed(get_faa_and_gff)(i, outDir_faa=info[1], outDir_gff=info[2]) for i in info[0])
        print(ide, ' done!')

    # for ide in ['galardini', 'botnar', 'colibafi', 'coliville', 'ecoref', 'lbc', 'roar', 'septicoli', 'touchon']:
    for ide in ['hong']:
        cmd = f'/nfs/home/smiravet/KTYPS_DEV/code/sh_utils/fetch_paths.sh /nfs/home/smiravet/KTYPS_DEV/scratch/raw/{ide}/{ide}_faa/ ".*\.faa$" /nfs/home/smiravet/KTYPS_DEV/scratch/raw/{ide}/fetch_annotations.txt'
        os.system(cmd)