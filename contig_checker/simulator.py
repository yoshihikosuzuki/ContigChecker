from random import choices
from BITS.utils import run_command


def generate_random_seq(length, bases=['A', 'C', 'G', 'T']):
    return ''.join(choices(bases, k=length))


def generate_tr_seq(tot_len, terminal_repeat_len):
    assert tot_len > 2 * terminal_repeat_len, "Invalid sequence lengths"
    tr_seq = generate_random_seq(terminal_repeat_len)
    u_seq = generate_random_seq(tot_len - 2 * terminal_repeat_len)
    return tr_seq + u_seq + tr_seq


def generate_reads(seq_fname, depth):
    prefix = f"depth{depth}"
    out_fname = f"{prefix}_0001.fastq"

    # Generate fastq using PBSIM
    command = ' '.join([f"pbsim",
                        f"--data-type CLR",
                        f"--depth {depth}",
                        f"--prefix {prefix}",
                        f"--sample-fastq ~/software/pbsim-1.0.3/sample/sample.fastq",   # TODO: parameterize
                        f"{seq_fname}"])
    run_command(command)

    # Convert fastq to daligner-acceptable fasta
    return run_command(f"gsed -n '2~4p' {out_fname}").strip().split('\n')


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Simulate a sequence with a pair of terminal repeats and noisy reads sampled from it using PBSIM")
    parser.add_argument('--seq_tot_len',
                        '-L',
                        type=int,
                        default=100000,
                        help="Overall sequence length to be generated.")
    parser.add_argument('--terminal_repeat_len',
                        '-l',
                        type=int,
                        default=2000,
                        help="Length of terminal repeats inside the sequence.")
    parser.add_argument('--read_depth',
                        '-d',
                        type=int,
                        default=30,
                        help="Read depth to be generaeted from the sequence.")
    args = parser.parse_args()

    true_seq = generate_tr_seq(args.seq_tot_len, args.terminal_repeat_len)
    with open("true_sequence.fasta", 'w') as f:
        f.write(f">Seq/0/0_{len(true_seq)}\n{true_seq}\n")

    reads = generate_reads("true_sequence.fasta", args.read_depth)
    with open("reads.fasta", 'w') as f:
        for i, read in enumerate(reads):
            f.write(f">Read/{i}/0_{len(read)}\n{read}\n")
