from collections import defaultdict, Counter
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import plotly.offline as py
import plotly.graph_objs as go

from BITS.utils import run_command

plt.style.use('ggplot')
py.init_notebook_mode()


def main(command_type="all",
         contigs_fasta=None,
         reads_fasta=None):
    pass


def load_args():
    import argparse
    parser = argparse.ArgumentParser(description="Call methylation motifs for each genome bin")
    parser.add_argument('command_type',
                        type=str,
                        help="Specify a module you want to run. {map, plot, categorize, modify}")
    parser.add_argument('contigs_fasta',
                        type=str,
                        default=None,
                        help="A fasta file of the contigs assembled.")
    parser.add_argument('reads_fasta',
                        type=str,
                        default=None,
                        help="A fasta file of the reads used for assembly.")
    return parser.parse_args()


if __name__ == "__main__":
    main(**load_args())
