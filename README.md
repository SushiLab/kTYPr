# 🧬 kTYPr

**kTYPr** is a command-line tool for predicting K-antigen type classifications in *Escherichia coli* using HMM-based annotation profiling. It supports both whole-genome and flanking-region prediction modes and can be run in parallel on multiple files.

---

## 📦 Installation

### Option 1 – Install via `pip` (recommended for most users)

If you already have Python and pip set up:

```bash
git clone https://github.com/SushiLab/kTYPr.git
cd kTYPr
pip install .
```

This will install `ktypr` and its dependencies.

To install in **editable/development mode** (so changes to the code apply immediately):

```bash
pip install -e .
```

Note this installation can be done within an existing conda (slower, can take minutes) or venv environment (faster, seconds).

### Option 2 – Use a Conda environment (safer, isolated setup)

If you use [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Anaconda](https://www.anaconda.com):

```bash
git clone https://github.com/SushiLab/kTYPr.git
cd kTYPr
conda env create -f environment.yml
conda activate ktypr_env
```

This will create and activate a new environment named `ktypr_env` with all dependencies and install `ktypr` in editable mode.

---

## 🚀 Usage

After installation, you can run the tool with:

```bash
ktypr --help
```

### Basic command

```bash
ktypr -i <input_path> 
```

### Required argument

* `-i`, `--input`:
  Path to a single file (e.g., `.faa`, `.fna`, `.gbk`), a directory containing multiple annotation files, or a `.txt` file listing full paths (one per line) to the genome/annotation files.

---

## 🔧 Options

| Option               | Description                                                                                                                                                                                                                                                                |
| -------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `-i`, `--input`      | **Input path**, which can be:<br>• A genome file (`.fasta`, `.gbk`, etc.)<br>• An annotation file (`.faa`)<br>• A directory containing such files<br>• A `.txt` file listing paths to input files. This is a required argument.                                            |
| `-o`, `--output`     | **Output directory** where all results will be saved. Defaults to `./ktypr_results`. The directory will be created if it does not exist.                                                                                                                                   |
| `-m`, `--mode`       | **Prediction mode**:<br>• `0` (default): *Flanking mode* — only genes located upstream and downstream of the **kpsC** gene, within a specified flanking window, are used for prediction.<br>• `1`: *Whole genome mode* — all annotated genes in the genome are considered. |
| `-f`, `--flank`      | **Flanking window size** in base pairs around the **kpsC** gene to evaluate (default: `30000`). This option is ignored when using whole genome mode (`-m 1`).                                                                                                              |
| `-n`, `--n-jobs`     | **Number of parallel jobs** to run:<br>• `-1` (default): use all available CPU cores.<br>• `1`: run sequentially without parallelism.<br>• Any other positive integer specifies the exact number of CPU cores to use.                                                      |
| `-p`, `--prefix`     | **Optional prefix** to prepend to all output file names. If not provided, the base name of each input genome or annotation file will be used as the prefix.                                                                                                                |
| `-s`, `--short`      | Flag to use metagenomic mode in prodigal gene calling in all input sequences. Recommended when short sequences are provided.                                                                                                                |
| `-r`, `--reannotate` | Flag to **force re-annotation** of genes using Prodigal, even if annotations are already present in the genome file. Useful to ensure consistent annotations when needed.         
| `-c`, `--clinker`    | Flag to produce [clinker](https://github.com/gamcil/clinker) reports. This can be computationally expensive, so it does not run by default.                                                                                                |
| `-v`, `--verbose`    | Enable **verbose mode** for detailed logging and debugging information during the run.                                                                                                                                                                                     |

### Genome annotation considerations

Please note the following aspects:
- kTYPr accepts annotated genomes (as `.faa` or `.gbk`), in this case, gene calling will not be run. Use `-r` to force the reannotation (only if genbank provided). 
- We use [pyrodigal](https://github.com/althonos/pyrodigal) in `single` mode for gene annotation by default. This can raise errors when the sequence is not long enough for training a gene caller specific for the input genome (min. length is 20,000 bases). In these cases, `meta` mode will be activated by default; take into account this can produce less accurate gene sequences (shorter/longer than they should). 
   - For a consistent gene calling, you can set `-s` to treat all provided sequences as short sequences and call genes with [pyrodigal](https://github.com/althonos/pyrodigal) metagenomic mode. 

---

## 🧪 Example

Run K-type prediction for whole-genomes (`--mode 1`) on a folder containing `.fa` files using all available cores (`-n -1`) and producing clinker reports (`-c`), saving results to a custom folder (`-o results/`) and printing in the terminal the running processes (`-v`):

```bash
ktypr -i ./test/genomes/fasta -o ./results/ --mode 1 -n -1 -c -v
```

Run in flanking mode using a text file of paths using a custom flanking size and with a custom prefix to name all the files:

```bash
ktypr -i genome_list.txt --flank 25000 -p ecoli_run
```

The expected runtime for a single genome prediction should be less than 30 seconds on a typical desktop computer.

---

## 🗃 Output

For each input genome, kTYPr creates in a folder per genome:

| Filename suffix         | Description                                                                                               |
| ----------------------- | --------------------------------------------------------------------------------------------------------- |
| `.gff`                  | Annotation file in GFF format. Contains coordinates and features of predicted genes.                      |
| `.faa`                  | Protein sequences of all annotated genes in FASTA format.                                                 |
| `_flanks.faa`           | Protein sequences extracted from the flanking region around the **kpsC** gene (only in flanking mode).    |
| `_hits.tsv.gz`          | Compressed TSV file listing all detected HMM hits (annotations) with their scores and locations.          |
| `_filtered_hits.tsv.gz` | Compressed TSV file with filtered HMM hits after applying score thresholds.                               |
| `_ktypr.tsv`            | Summary TSV file containing the final K-antigen type prediction results for the genome or annotation set. |
| `.gbk`                  | Full genome file with annotations in GenBank format, optionally including re-annotation results.          |
| `clink.html`            | [Clinker](https://github.com/gamcil/clinker) HTML report against the best K-antigen type predicted.        |

### Collection results

A summary table of all classifications is also created if more than one genome is given as `<prefix>results_ktypr.tsv` in the selected output directory. This file provides a summary of the K-antigen (K-type) prediction results for each genome or annotation set. It contains both the final predicted type and detailed match statistics for conserved regions and all candidate K-types.

Columns description:

| Column name            | Description                                                                                                           |
| ---------------------- | --------------------------------------------------------------------------------------------------------------------- |
| `predicted`            | **Best-matching K-type** predicted for the genome based on gene content and HMM match scores.                         |
| `pred_nr_genes`        | Number of genes expected for the predicted K-type (based on its reference profile).                                   |
| `pred_genes_in_genome` | Number of those expected genes actually found in the genome.                                                          |
| `pred_acc_bitscore`    | Sum of HMM bitscores across all hits matching the predicted K-type. A proxy for overall match strength.               |
| `pred_is_complete`     | Indicates if all expected genes for the predicted K-type were found in the genome (`1` = complete, `0` = incomplete). |

These are followed by conserved K-locus genes shared across many K-types, such as those in the KpsEDCSMT and KpsFU operons:

| Column name                 | Description                                                     |
| --------------------------- | --------------------------------------------------------------- |
| `KpsEDCSMT_nr_genes`        | Number of conserved EDCSMT genes expected.                      |
| `KpsEDCSMT_genes_in_genome` | Number of EDCSMT genes found in the genome.                     |
| `KpsEDCSMT_acc_bitscore`    | Cumulative bitscore for the EDCSMT gene matches.                |
| `KpsEDCSMT_is_complete`     | Whether all EDCSMT genes were found (`1`) or not (`0`).         |
| `KpsFU_nr_genes`, ...       | Same structure as above, but for the conserved **FU** gene set. |

For further exploration, and forF each known K-type in our database  the following columns provide detailed information:

| Column suffix             | Description                                                           |
| ------------------------- | --------------------------------------------------------------------- |
| `<KTYPE>_nr_genes`        | Number of genes expected for this K-type.                             |
| `<KTYPE>_genes_in_genome` | Number of expected genes found in the genome.                         |
| `<KTYPE>_acc_bitscore`    | Sum of HMM match bitscores for this K-type’s genes.                   |
| `<KTYPE>_is_complete`     | Whether the K-type is fully present in the genome (`1`) or not (`0`). |

As example of the collection output, you can find in [test/output](test/output) the results of running:

```bash
ktypr -i ./test/genomes/fasta -v -o ./test/output/ -p run1_
```

---

## 📖 Citation
If you use kTYPr in your research, please cite:

Roese. et al. XXXXX. Journal (2025). LINK

---

## ❓ Need Help?

```bash
ktypr --help
```

or contact the maintainer at *[smiravet@ethz.ch](mailto:smiravet@ethz.ch)*
