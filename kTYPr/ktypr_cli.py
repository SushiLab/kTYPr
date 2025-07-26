import argparse
from pathlib import Path

def parse_args():
    parser = argparse.ArgumentParser(description="Run kTYPr on a genome or annotation set.")

    parser.add_argument("-i", "--input", required=True,
                        help="Input genome file (fasta/genbank), annotation file (fna/faa), directory of such files, or a .txt file listing paths"
    )

    parser.add_argument("-o", "--output", default="./ktypr_results",
                        help="Output directory (default: ./ktypr_results)"
    )

    parser.add_argument("-m", "--mode", type=int, choices=[0, 1], default=0,
                        help=(
                            "Prediction mode:\n"
                            "  0 = flanking (default): use genes upstream and downstream of the kpsC gene,\n"
                            "      within a window defined by the --flank parameter (in base pairs)\n"
                            "  1 = whole genome: use all annotated genes across the genome"
                        )
    )

    parser.add_argument("-f", "--flank", type=int, default=30000,
                        help="Flanking window size in base pairs to extract and only evaluate genes around kpsC (default: 30000). This is not considered when using whole-genome mode."
    )

    parser.add_argument("-n", "--n-jobs", type=int, default=-1,
                        help="Number of parallel jobs (default: -1 = use all cores, 1 = disable parallelism, any other number = number of cores to use)"
    )

    parser.add_argument("-p", "--prefix", default=None, 
                        help="Optional prefix to add to output file names. By default, use the genome/annotation file basename"
    )

    parser.add_argument("-r", "--reannotate", action="store_true",
                        help="Reannotate genes using Prodigal, even if annotations are present in the GenBank file."
    )

    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Enable verbose output for debugging purposes."
    )
    
    return parser.parse_args()

def main():
    args = parse_args()
    input_path = Path(args.input)
    output_dir = Path(args.output)
    
    # Determine if parallelism should be enabled
    parallel = args.n_jobs != 1

    # Call internal function (you must define this accordingly)
    run_ktypr(
        input_path=input_path,
        output_dir=output_dir,
        mode=args.mode,
        flank=args.flank,
        parallel=parallel,
        n_jobs=args.n_jobs,
        prefix=args.prefix
    )

if __name__ == "__main__":
    main()
