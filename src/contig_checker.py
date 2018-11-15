from logzero import logger
from .cc_map import run_mapping
from .cc_annotate import annotate_regions
from .cc_plot import plot_figure
from .cc_modify import modify_sequence


def main():
    args = load_args()
    if args.command_type == "map":
        run_mapping(args.contigs_fasta, args.reads_fasta, args.n_core)
    elif args.command_type == "annotate":
        annotate_regions()
    elif args.command_type == "plot":
        plot_figure()
    elif args.command_type == "modify":
        modify_sequence()
    else:
        logger.error(f"No module '{args.command_type}' is offered.")


def load_args():
    import argparse
    parser = argparse.ArgumentParser(description="Call methylation motifs for each genome bin")
    parser.add_argument('command_type',
                        type=str,
                        help="Specify a module you want to run. {map, plot, categorize, modify}")
    parser.add_argument('--contigs_fasta',
                        '-c',
                        type=str,
                        default=None,
                        help="A fasta file of the contigs assembled.")
    parser.add_argument('--reads_fasta',
                        '-r',
                        type=str,
                        default=None,
                        help="A fasta file of the reads used for assembly.")
    parser.add_argument('--n_core',
                        '-n',
                        type=int,
                        default=1,
                        help="Degree of job parallelization. [1]")
    parser.add_argument('--plot_read_id',
                        '-i',
                        type=str,
                        default=None,
                        help="ID of a read to be plotted in the 'plot' mode.")
    return parser.parse_args()


if __name__ == "__main__":
    main()
