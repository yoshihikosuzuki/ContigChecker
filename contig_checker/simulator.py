import random
from BITS.utils import run_command


def generate_random_seq(length, bases=['A', 'C', 'G', 'T']):
    return ''.join(random.choices(bases, k=length))


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


def sample_seq(s, mu=8.9, sigma=0.5):
    """
    Cut out a substring from <s> whose length obeys log-normal distribution.
    This function generates terminal-including sequences when they reach it during sampling.
    """
    direction = random.choice((-1, 1))   # NOTE: this is not strand, but for generating terminal-ending reads
    start = random.randint(0, len(s) - 1)
    length = int(random.lognormvariate(mu, sigma))
    return s[start:min(start + length, len(s))] if direction == 1 else s[max(0, start - length + 1):start + 1]


def introduce_noise(s):
    pass


def single_to_multi(s, width=100):
    """
    Cut a single sequence at every <width> bp and return as a list.
    """
    return [s[i:i + width] for i in range(0, len(s), width)]


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
        f.write(f">Seq/0/0_{len(true_seq)}\n")
        seqs = single_to_multi(true_seq)
        for seq in seqs:
            f.write(f"{seq}\n")

    """
    reads = generate_reads("true_sequence.fasta", args.read_depth)
    with open("reads.fasta", 'w') as f:
        for i, read in enumerate(reads):
            f.write(f">Read/{i}/0_{len(read)}\n")
            seqs = single_to_multi(read)
            for seq in seqs:
                f.write(f"{seq}\n")
    """

    with open("reads.fasta", 'w') as f:
        length_sum = 0
        index = 0
        while length_sum / len(true_seq) < args.read_depth:
            read = sample_seq(true_seq)
            f.write(f">Read/{index}/0_{len(read)}\n")   # TODO: rewrite using BITS
            seqs = single_to_multi(read)
            for seq in seqs:
                f.write(f"{seq}\n")
            index += 1
            length_sum += len(read)
