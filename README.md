# 🧬 kTYPr

**kTYPr** is a command-line tool for predicting K-antigen type classifications in bacterial genomes using HMM-based annotation profiling. It supports both whole-genome and flanking-region prediction modes and can be run in parallel on multiple files.

---

## 📦 Installation

You can install the tool by cloning the repository and using `pip`:

```bash
git clone https://github.com/SushiLab/kTYPr.git
cd kTYPr
pip install .
```

After installation, you can run the tool with:

```bash
ktypr --help
```

---

## 🚀 Usage

### Basic command

```bash
ktypr -i <input_path> [options]
```

### Required argument

* `-i`, `--input`:
  Path to a single file (e.g., `.faa`, `.fna`, `.gbk`), a directory containing multiple annotation files, or a `.txt` file listing full paths (one per line) to the genome/annotation files.

---

## 🔧 Options

| Option           | Description                                                                                                                                                                                   |
| ---------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `-o`, `--output` | Output directory (default: `./ktypr_results`)                                                                                                                                                 |
| `-m`, `--mode`   | Prediction mode:<br> `0`: Flanking (default) – only genes around the **kpsC** gene within a given window are considered.<br> `1`: Whole genome – all annotated genes are used for prediction. |
| `-f`, `--flank`  | Flanking region size in base pairs to evaluate around **kpsC** (default: `30000`). Ignored in whole-genome mode.                                                                              |
| `-n`, `--n-jobs` | Number of parallel jobs (default: `-1` = all cores, `1` = sequential)                                                                                                                         |
| `-p`, `--prefix` | Optional prefix to add to the output filenames. If not set, the base name of each genome file is used.                                                                                        |

---

## 🧪 Example

Run K-type prediction on a folder of `.faa` files using all available cores, saving results to a custom folder:

```bash
ktypr -i data/genomes_faa/ -o results/ --mode 1 -n -1
```

Run in flanking mode using a text file of paths:

```bash
ktypr -i data/genomes_list.txt --mode 0 --flank 25000 -p ecoli_run
```

---

## 🗃 Output

For each input genome, kTYPr creates:

* A gzipped TSV file of HMM hits (`*_hits.tsv.gz`)
* A TSV file with the final K-type classification (`*_ktyps.tsv`)

A summary table of all classifications is also created if a list or directory is given.

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