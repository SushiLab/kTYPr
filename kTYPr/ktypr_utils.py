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
        full_genome_seq = b''.join([bytes(record.seq) for record in records])
        orf_finder = pyrodigal.GeneFinder()
        orf_finder.train(full_genome_seq)
    return orf_finder

def get_faa_and_gff(fil, outDir_faa, outDir_gff, meta=False):
    """ Runs prodigal on a genome file. It accepts both compressed (as gz) and uncompressed fasta files """

    ide = fil.split('/')[-1].split('.')[0]
    faa_basename = f'{outDir_faa}/{ide}.faa'
    gff_basename = f'{outDir_gff}/{ide}.gff'
    ffo = open(faa_basename, "w")
    gfo = open(gff_basename, "w")

    if fil.endswith('.gz'):
        with gzip.open(fil, "rt") as f:
            records = list(SeqIO.parse(f, "fasta"))  # <-- List once
    else:
        with open(fil, "rt") as f:
            records = list(SeqIO.parse(f, "fasta"))

    print(f'Processing {ide} with {len(records)} records')

    orf_finder = set_pyrodigal_model(records, meta=meta)    
    for i, record in enumerate(records):  # reuse list
        print(i, ide)
        genes = orf_finder.find_genes(bytes(record.seq))
        genes.write_gff(gfo, sequence_id=f'{ide}__{record.id}', header=(i == 0))
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


def set_pyrodigal_model(records, meta=False):
    """ Sets up the Pyrodigal gene finder model based on the provided records."""
    if meta:
        print("Setting up Pyrodigal in metagenomic mode")
        return pyrodigal.GeneFinder(meta=True)
    else:
        print("Setting up Pyrodigal in genomic mode")
        full_genome_seq = b''.join([bytes(record.seq) for record in records])
        orf_finder = pyrodigal.GeneFinder()
        orf_finder.train(full_genome_seq)
        return orf_finder

def split_multigenbank_to_faa_and_gff(input_genbank, outDir_faa, outDir_gff, reannotate=False, meta=False):
    """
    Splits a multi-GenBank file into individual genome files and extracts:
    - A `.faa` file with the CDS (coding sequences) in protein format.
    - A `.gff` file with annotations for the genes.

    If reannotate is True or if no CDS are found in a record, Pyrodigal will be used to predict genes.

    Args:
        input_genbank (str): Path to the multi-GenBank file.
        outDir_faa (str): Path to the directory to store .faa files.
        outDir_gff (str): Path to the directory to store .gff files.
        reannotate (bool): Whether to reannotate with Pyrodigal.
        meta (bool): Whether to use metagenomic mode for Pyrodigal.
    """
    outDir_faa = Path(outDir_faa)
    outDir_gff = Path(outDir_gff)
    outDir_faa.mkdir(parents=True, exist_ok=True)
    outDir_gff.mkdir(parents=True, exist_ok=True)

    records = list(SeqIO.parse(input_genbank, "genbank"))

    # Set up Pyrodigal only once if needed
    orf_finder = None
    if reannotate or any(not any(f.type == "CDS" for f in rec.features) for rec in records):
        orf_finder = set_pyrodigal_model(records, meta=meta)
    
    for i, record in enumerate(records):
        genome_id = record.name
        faa_file = outDir_faa / f"{genome_id}.faa"
        gff_file = outDir_gff / f"{genome_id}.gff"

        with open(faa_file, "w") as faa_out, open(gff_file, "w") as gff_out:
            # Write GFF header
            gff_out.write("##gff-version 3\n")
            gff_out.write(f"##sequence-region {genome_id} 1 {len(record)}\n")

            cds_features = [f for f in record.features if f.type == "CDS"]

            if reannotate or not cds_features:
                # Run Pyrodigal gene prediction
                genes = orf_finder.find_genes(bytes(record.seq))
                genes.write_gff(gff_out, sequence_id=genome_id, header=False)
                genes.write_translations(faa_out, sequence_id=genome_id)
            else:
                for feature in cds_features:
                    cds_seq = feature.extract(record.seq)
                    protein_seq = cds_seq.translate(to_stop=True)

                    protein_id = feature.qualifiers.get("locus_tag", ["unknown"])[0]
                    gene_id = feature.qualifiers.get("gene", ["unknown"])[0]

                    faa_record = SeqRecord(protein_seq, id=gene_id + '_' + protein_id, description="")
                    SeqIO.write(faa_record, faa_out, "fasta")

                    start = int(feature.location.start) + 1
                    end = int(feature.location.end)
                    strand = "+" if feature.location.strand == 1 else "-"
                    attributes = f"ID={protein_id};"
                    if "product" in feature.qualifiers:
                        attributes += f"product={feature.qualifiers['product'][0]};"
                    gff_line = f"{genome_id}\tGenBank\tCDS\t{start}\t{end}\t.\t{strand}\t0\t{attributes}\n"
                    gff_out.write(gff_line)
                    
# I/O

def process_genome_file(genome_path, outDir_faa, outDir_gff, reannotate=False, meta=False):
    genome_path = Path(genome_path)
    if genome_path.suffix in {".gb", ".gbk", ".gbff"}:
        # GenBank format
        split_multigenbank_to_faa_and_gff(str(genome_path), outDir_faa, outDir_gff, reannotate=reannotate)
    elif genome_path.suffix in {".fna", ".fa", ".fasta"}:
        # FASTA format
        get_faa_and_gff(str(genome_path), outDir_faa=outDir_faa, outDir_gff=outDir_gff, meta=meta)
    else:
        print(f"Unknown file format for {genome_path}, skipping.")

def prepare_input(inFile, outDir, prefix='', extract_annotations=True, reannotate=False, n_jobs=30, meta=False):
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
        # List them in the same file
        annotations_list = os.path.join(outDir, "fetch_annotations.txt")
        with open(annotations_list, "w") as f:
            for path in input_paths:
                f.write(path + "\n")

    elif extract_annotations:
        basedir = os.path.dirname(os.path.abspath(inFile))
        ide = os.path.basename(basedir)

        # fetch_txt_path = os.path.join(basedir, "fetch.txt")
        # with open(fetch_txt_path, "w") as f:
        #     for path in input_paths:
        #         f.write(path + "\n")

        faa_o = os.path.join(outDir, f"{prefix}{ide}_faa")
        gff_o = os.path.join(outDir, f"{prefix}{ide}_gff")
        os.makedirs(faa_o, exist_ok=True)
        os.makedirs(gff_o, exist_ok=True)

        # Run gene calling
        Parallel(n_jobs=n_jobs)(
            delayed(process_genome_file)(genome_path, outDir_faa=faa_o, outDir_gff=gff_o, reannotate=reannotate, meta=meta)
            for genome_path in input_paths
        )

        annotations_list = os.path.join(outDir, "fetch_annotations.txt")
        cmd = f'{SCRIPT_DIR}/fetch_paths.sh {faa_o} ".*\\.faa$" {annotations_list}'
        os.system(cmd)
    else:
        raise ValueError("Input does not appear to be annotations, and annotation extraction is disabled.")
    
    return annotations_list

#import os
#import glob
#from joblib import Parallel, delayed
#from pathlib import Path
#
#SCRIPT_DIR = Path(__file__).resolve().parent  # Adjust this if SCRIPT_DIR is defined elsewhere
#
#def prepare_input(inFile, outDir, prefix='', extract_annotations=True, reannotate=False, n_jobs=30):
#    """
#    Prepares a list of annotation files from the given input:
#    - If input is annotations (.faa), return directly.
#    - If input is genome files (FASTA/GenBank), extract annotations (optionally using a user-defined prefix).
#    
#    Returns:
#        annotations_list (str): Path to file listing all .faa annotation paths.
#        tracking_dict (dict): Dictionary with keys 'ides', 'path_faa', 'path_gff'
#    """
#    def resolve_paths(inFile):
#        if inFile.endswith(".txt"):
#            with open(inFile) as f:
#                return [line.strip() for line in f if line.strip()]
#        elif os.path.isdir(inFile):
#            extensions = ("*.fasta", "*.fa", "*.fna", "*.gb", "*.gbk", "*.faa")
#            files = []
#            for ext in extensions:
#                files.extend(glob.glob(os.path.join(inFile, ext)))
#            return sorted(files)
#        elif os.path.isfile(inFile):
#            return [inFile]
#        else:
#            raise ValueError(f"Invalid input: {inFile} is neither a file nor a directory.")
#
#    input_paths = resolve_paths(inFile)
#    tracking_dict = {'ides': [], 'path_faa': [], 'path_gff': []}
#
#    if all(path.endswith(".faa") for path in input_paths):
#        # Case: annotations already exist
#        annotations_list = os.path.join(outDir, "fetch_annotations.txt")
#        with open(annotations_list, "w") as f:
#            for path in input_paths:
#                f.write(path + "\n")
#                ide = os.path.splitext(os.path.basename(path))[0]
#                tracking_dict['ides'].append(ide)
#                tracking_dict['path_faa'].append(path)
#                tracking_dict['path_gff'].append(None)
#
#    elif extract_annotations:
#        basedir = os.path.dirname(os.path.abspath(inFile))
#        ide = os.path.basename(basedir)
#
#        faa_o = os.path.join(outDir, f"{prefix}{ide}_faa")
#        gff_o = os.path.join(outDir, f"{prefix}{ide}_gff")
#        os.makedirs(faa_o, exist_ok=True)
#        os.makedirs(gff_o, exist_ok=True)
#
#        # Run gene calling
#        Parallel(n_jobs=n_jobs)(
#            delayed(process_genome_file)(
#                genome_path,
#                outDir_faa=faa_o,
#                outDir_gff=gff_o,
#                reannotate=reannotate
#            )
#            for genome_path in input_paths
#        )
#
#        # Track outputs
#        for genome_path in input_paths:
#            base = os.path.splitext(os.path.basename(genome_path))[0]
#            faa_file = os.path.join(faa_o, base + ".faa")
#            gff_file = os.path.join(gff_o, base + ".gff")
#            tracking_dict['ides'].append(base)
#            tracking_dict['path_faa'].append(faa_file)
#            tracking_dict['path_gff'].append(gff_file)
#
#        annotations_list = os.path.join(outDir, "fetch_annotations.txt")
#        cmd = f'{SCRIPT_DIR}/fetch_paths.sh {faa_o} ".*\\.faa$" {annotations_list}'
#        os.system(cmd)
#
#    else:
#        raise ValueError("Input does not appear to be annotations, and annotation extraction is disabled.")
#    
#    return annotations_list, tracking_dict
#
#


from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
import gzip

def parse_gff(gff_file):
    """
    Minimal GFF3 parser to extract features.
    Yields tuples: (seq_id, source, type, start, end, strand, attributes_dict)
    """
    open_func = gzip.open if gff_file.endswith(".gz") else open
    with open_func(gff_file, "rt") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) != 9:
                continue
            seqid, source, type_, start, end, score, strand, phase, attributes = parts
            attr_dict = {}
            for attr in attributes.split(";"):
                if "=" in attr:
                    key, value = attr.split("=", 1)
                    attr_dict[key] = value
            yield seqid, source, type_, int(start), int(end), strand, attr_dict


def create_genbank_from_inputs(input_file, gff_file, genome_fasta, selected_ids, output_genbank, id_attribute="ID", from_annotations=False):
    """
    Create a GenBank file with selected genes from either:
    - annotations (.faa): builds CDS-only GenBank
    - genome input (.fa/.gb + .gff): builds full GenBank with selected features
    """
    output_genbank = Path(output_genbank)

    if from_annotations:
        # From .faa input → CDS-only GenBank
        records = list(SeqIO.parse(input_file, "fasta"))
        filtered_records = [r for r in records if any(gid in r.id for gid in selected_ids)]
        for r in filtered_records:
            r.annotations["molecule_type"] = "DNA"
        SeqIO.write(filtered_records, output_genbank, "genbank")

    else:
        # Parse genome fasta
        genome = SeqIO.to_dict(SeqIO.parse(genome_fasta, "fasta"))

        # Build new SeqRecords with selected features
        record_dict = {k: SeqRecord(seq=v.seq, id=v.id, name=v.name, description="", annotations={"molecule_type": "DNA"})
                       for k, v in genome.items()}

        for seqid, source, type_, start, end, strand, attr in parse_gff(gff_file):
            if id_attribute not in attr:
                continue
            feature_id = attr[id_attribute]
            if isinstance(feature_id, str):
                match = any(gid == feature_id for gid in selected_ids)
            else:
                match = any(gid in feature_id for gid in selected_ids)

            if not match or seqid not in record_dict:
                continue

            location = FeatureLocation(start - 1, end, strand=1 if strand == "+" else -1)
            qualifiers = {k: [v] for k, v in attr.items()}
            feat = SeqFeature(location=location, type=type_, qualifiers=qualifiers)
            record_dict[seqid].features.append(feat)

        # Write only records with features
        filtered_records = [r for r in record_dict.values() if r.features]
        SeqIO.write(filtered_records, output_genbank, "genbank")


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