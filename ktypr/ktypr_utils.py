import os
import glob
import gzip
import shutil
import zipfile
import pyrodigal
import pandas as pd
from . import ktypr_hmms as kh
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from pathlib import Path
from collections import defaultdict


### GLOBAL VARIABLES
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

### Gene calling from faa and gbk functions

def set_pyrodigal_model(records, meta=False, verbose=True):
    """ Sets up the Pyrodigal gene finder model based on the provided records."""
    if meta:
        if verbose:
            print("Setting up Pyrodigal in metagenomic mode")
        return pyrodigal.GeneFinder(meta=True)
    else:
        if verbose:
            print("Setting up Pyrodigal in genomic mode")
        full_genome_seq = b''.join([bytes(record.seq) for record in records])
        if len(full_genome_seq) < 20000:
            return set_pyrodigal_model(records=None, meta=True, verbose=verbose)           # To avoid <20kb error
        else:
            orf_finder = pyrodigal.GeneFinder()
            orf_finder.train(full_genome_seq)
            return orf_finder


def get_faa_and_gff(fil, outfaa, outgff, ide, meta=False, verbose=True):
    """
    Runs Prodigal on a genome file (compressed or not).
    Writes translated proteins to `outfaa` and gene annotations to `outgff`.
    """
    open_func = gzip.open if str(fil).endswith('.gz') else open
    with open_func(fil, "rt") as f, open(outfaa, "w") as ffo, open(outgff, "w") as gfo:
        records = list(SeqIO.parse(f, "fasta"))
        orf_finder = set_pyrodigal_model(records, meta=meta, verbose=verbose)
        if verbose:
            print(f'Gene calling running on {ide} with {len(records)} records')
        for i, record in enumerate(records):
            seq_id = f"{ide}__{record.id}"
            genes = orf_finder.find_genes(bytes(record.seq))
            genes.write_gff(gfo, sequence_id=seq_id, header=(i == 0))
            genes.write_translations(ffo, sequence_id=seq_id)


def split_genbank_to_faa_and_gff(input_genbank, outfaa, outgff, ide, reannotate=False, meta=False, verbose=True):
    """
    Splits a multi-GenBank file into individual genome files and extracts:
    - A `.faa` file with the CDS (coding sequences) in protein format.
    - A `.gff` file with annotations for the genes.

    If reannotate is True or if no CDS are found in a record, Pyrodigal will be used to predict genes.

    Args:
        input_genbank (str): Path to the multi-GenBank file.
        outfaa (str): Path to the .faa file.
        outgff (str): Path to the .gff file.
        reannotate (bool): Whether to reannotate with Pyrodigal.
        meta (bool): Whether to use metagenomic mode for Pyrodigal.
    """
    records = list(SeqIO.parse(input_genbank, "genbank"))

    # Set up Pyrodigal only once if needed
    orf_finder = None
    if reannotate or any(not any(f.type == "CDS" for f in rec.features) for rec in records):
        orf_finder = set_pyrodigal_model(records, meta=meta)

    # Check number of records
    flag = False
    if len(records) == 1:
        flag = True
        record = records[0]
    else:
        if verbose:
            if len(records) == 0:
                print(f"No records found in {input_genbank}. Skipping.")
            else:
                print(f"Multiple records ({len(records)}) in {input_genbank}, please provide a single record GenBank file.")
    
    # Split records to faa and gff files
    if flag:
        with open(outfaa, "w") as faa_out, open(outgff, "w") as gff_out:
            # Write GFF header
            gff_out.write("##gff-version 3\n")
            gff_out.write(f"##sequence-region {ide} 1 {len(record)}\n")
            cds_features = [f for f in record.features if f.type == "CDS"]
            if reannotate or not cds_features:
                # Run Pyrodigal gene prediction
                genes = orf_finder.find_genes(bytes(record.seq))
                genes.write_gff(gff_out, sequence_id=ide, header=False)
                genes.write_translations(faa_out, sequence_id=ide)
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
                    gff_line = f"{ide}\tGenBank\tCDS\t{start}\t{end}\t.\t{strand}\t0\t{attributes}\n"
                    gff_out.write(gff_line)

# Parsers and gene subsetting

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


# Processing input 

def check_valid_input(files, extensions, verbose=True):
    """
    Checks if files exist and match the given extensions.

    Args:
        files (list): List of file paths.
        extensions (list): List of allowed extensions (e.g., ['.fa', '.gbk']).
        verbose (bool): Print messages for invalid files.

    Returns:
        list: Valid file paths.
    """
    # Normalize extensions to start with dot and lowercase
    norm_exts = {ext if ext.startswith('.') else f'.{ext}' for ext in extensions}
    norm_exts = {ext.lower() for ext in norm_exts}

    valid = []
    invalid = []
    for file in files:
        p = Path(file)
        if p.is_file() and p.suffix.lower() in norm_exts:
            valid.append(str(p))
        else:
            invalid.append(str(p))
    if verbose:
        for file in invalid:
            print(f"Invalid file: {file}. Expected one of {sorted(norm_exts)}.")
        if not valid:
            print("No valid files found.")
    return valid


def resolve_paths(inFile, extensions = [".fasta", ".fa", ".fna", ".gb", ".gbk", ".faa"], verbose=True):
    """
    Resolves input paths based on the type of input provided.
    - If a .txt file, reads paths from it.
    - If a directory, collects all fasta/genbank files.
    - If a single file, returns that file.
    """

    # Get and check valid if list of files as txt

    if inFile.endswith(".txt"):
        with open(inFile) as f:
            files = [line.strip() for line in f if line.strip()]
    elif os.path.isdir(inFile):
        files = []
        for ext in extensions:
            files.extend(glob.glob(os.path.join(inFile, f'*{ext}')))
    elif os.path.isfile(inFile):
        files = [inFile]
    else:
        raise ValueError(f"Invalid input: {inFile} is neither a file nor a directory.")
    
    return sorted(check_valid_input(files, extensions, verbose=verbose))


def process_genome_file(genome_path, outfaa, outgff, ide, reannotate=False, meta=False, verbose=True):
    genome_path = Path(genome_path)
    if genome_path.suffix in {".gb", ".gbk", ".gbff"}:
        # GenBank format
        split_genbank_to_faa_and_gff(str(genome_path), outfaa, outgff, ide, reannotate=reannotate, meta=meta, verbose=verbose)
    elif genome_path.suffix in {".fna", ".fa", ".fasta"}:
        # FASTA format
        get_faa_and_gff(str(genome_path), outfaa, outgff, ide, meta=meta, verbose=verbose)
    else:
        if verbose:
            print(f"Unknown file format for {genome_path}, skipping.")


def prepare_single_input(genome_path, outDir, prefix='', 
                         extract_annotations=True, reannotate=False, meta=False, verbose=True):
    """
    Checks and annotates with pyrodigal (if needed/selected) a single genome or annotation file.
    Returns a dictionary to handle results with identifiers, paths to .faa and .gff files,
    and expected outputs for hits and classification.
    """
    genome_path = Path(genome_path)
    identifier  = prefix + genome_path.stem
    outdir = Path(outDir) / identifier
    outdir.mkdir(parents=True, exist_ok=True)

    result = {
        'ide': identifier,
        'genome_path': str(genome_path),
        'outdir': str(outdir),
        'faa_path': str(outdir / f"{identifier}.faa"),
        'faa_flank_path': str(outdir / f"{identifier}_flanks.faa"),
        'gff_path': str(outdir / f"{identifier}.gff"),
        'gbk_path': str(outdir / f"{identifier}.gbk"),
        'hits': str(outdir / f"{identifier}_hits.tsv.gz"),
        'filt_hits': str(outdir / f"{identifier}_filtered_hits.tsv.gz"),
        'classification': str(outdir / f"{identifier}_ktypr.tsv"),
        'clinker': str(outdir / f"{identifier}_clink.html"),
        'annotated': 0,
        'done': 0
    }

    # Already annotated (input is a .faa file)
    if genome_path.suffix == ".faa":
        if verbose:
            print(f"A set of proteins (.faa) has been provided: {genome_path}")
        result['faa_path']  = str(genome_path)
        result['gff_path']  = None
        result['gbk_path']  = None
        result['clinker']   = None
        result['annotated'] = 1
        return result

    if not extract_annotations:
        raise ValueError(f"Annotation extraction disabled, but input is not .faa: {genome_path}")

    # Run annotation
    process_genome_file(
        genome_path,
        outfaa=result['faa_path'],
        outgff=result['gff_path'],
        ide=identifier,
        reannotate=reannotate,
        meta=meta,
        verbose=verbose
    )
    result['annotated'] = 1

    return result

# Processing output 

def parse_gff(gff_file):
    """
    Minimal GFF3 parser to extract features.
    Yields: (seqid, source, type_, start, end, strand, attributes_dict)
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


def create_genbank_from_inputs(result_dict, id_attribute="ID", from_annotations=False, verbose=True):
    """
    Create a GenBank file with selected genes from either:
    - From annotations (.faa): builds CDS-only GenBank
    - From genome input (.fa/.gb + .gff): builds full GenBank with selected features, including /translation
    """
    # Load translations from .faa
    translation_dict = {
        rec.id: str(rec.seq)
        for rec in SeqIO.parse(result_dict['faa_path'], "fasta")
    }

    # Load hits file and build feature metadata dictionary
    if result_dict['max_hits'] is not None and not result_dict['max_hits'].empty:
        hits_df  = result_dict['max_hits']
        hits_dict = (
            hits_df
            .groupby('query')['subject']
            .agg(lambda x: ';'.join(map(str, x)))
            .to_dict()
            )    
        selected_ids = set(hits_dict.keys())

        if from_annotations:
            # Filter only selected translations
            filtered_records = [
                rec for rec in SeqIO.parse(result_dict['faa_path'], "fasta")
                if any(gid in rec.id for gid in selected_ids)
            ]
            for rec in filtered_records:
                rec.annotations["molecule_type"] = "DNA"
            if not filtered_records and verbose:
                print("Warning: No matching records found in .faa input.")
            SeqIO.write(filtered_records, result_dict['gbk_path'], "genbank")

        else:
            # Determine file format and parse
            genome_path = result_dict['genome_path']
            genome_suffix = Path(genome_path).suffix.lower()
            fmt = "genbank" if genome_suffix in [".gbk", ".gb", ".genbank"] else "fasta"
            genome = SeqIO.to_dict(SeqIO.parse(genome_path, fmt))

            record_dict = {
                f"{result_dict['ide']}__{k}": SeqRecord(
                    seq=v.seq,
                    id=v.id,
                    name=v.name,
                    description="",
                    annotations={"molecule_type": "DNA"}
                )
                for k, v in genome.items()
            }

            # Parse GFF and build features
            for seqid, source, type_, start, end, strand, attr in parse_gff(result_dict['gff_path']):
                if id_attribute not in attr:
                    continue
                feature_id = attr[id_attribute]
                if not isinstance(feature_id, str):
                    continue
                if feature_id not in selected_ids or seqid not in record_dict:
                    continue

                location = FeatureLocation(start - 1, end, strand=1 if strand == "+" else -1)
                qualifiers = {k: [v] for k, v in attr.items()}

                # Add /translation if available
                if feature_id in translation_dict:
                    qualifiers["translation"] = [translation_dict[feature_id]]

                # Add /gene, /locus_tag, /protein_id if available
                if feature_id in hits_dict:
                    for key in ["gene", "locus_tag", "protein_id"]:
                        qualifiers[key] = [hits_dict[feature_id]]

                feature = SeqFeature(location=location, type=type_, qualifiers=qualifiers)
                record_dict[seqid].features.append(feature)

            # Write only records that have features
            filtered_records = [r for r in record_dict.values() if r.features]
            if not filtered_records and verbose:
                print("Warning: No matching features found. GenBank file will be empty.")
            SeqIO.write(filtered_records, result_dict['gbk_path'], "genbank")
    else:
        if verbose:
            print(f"No genbank produced for {result_dict['genome_path']}. No KpsC found, try in whole-genome mode.")


# CLINKER 

def get_clinker(results, verbose=True):
    """
    Extract specific K-antigen clusters from a zipped archive of GenBanks,
    which can then be used with clinker.

    `results` can be:
      - a dict with keys 'best' and 'clinker_dir'
      - a list of such dicts

    The function looks for GenBank files in the archive whose name contains
    the `best` identifier, and extracts them into the corresponding `clinker_dir`.
    """
    if isinstance(results, dict):
        results = [results]  # make it a list of one item

    # Extract sequentially
    archive_path = f"{SCRIPT_DIR}/data/reference_clusters.zip"
    with zipfile.ZipFile(archive_path, 'r') as zipf:
        for res in results:
            ks = res['ktype']
            out_dir = res['outdir']

            # Ensure the output directory exists
            os.makedirs(out_dir, exist_ok=True)

            # Find matching files and extract
            matched = [file for file in zipf.namelist() if f'{ks}_' in file]
            for file in matched:
                zipf.extract(file, path=out_dir)
                if verbose:
                    print(f"Extracted: {file} to {out_dir}")
    
    # Run the html production with clinker
    for res in results:
        original_gbk  = res['gbk_path']
        reference_dir = Path(res['outdir']) / "reference_clusters"
        if res['gff_path']:
            if verbose:
                print(f'Running clinker on {original_gbk}')

            # Copy and run clinker
            shutil.copy(original_gbk, reference_dir)
            cmd = f'clinker {reference_dir}/*.gbk -p {res["clinker"]}'

            if verbose:
                os.system(cmd)
                print(f'Clinker report in {res["clinker"]}')
            else:
                os.system(f"{cmd} > /dev/null 2>&1")

            # Safely remove 
            try:
                if reference_dir.exists() and reference_dir.is_dir():
                    shutil.rmtree(reference_dir)
                    if verbose:
                        print(f"Removed folder: {reference_dir}")
            except Exception as e:
                print(f"Failed to remove folder {reference_dir}: {e}")