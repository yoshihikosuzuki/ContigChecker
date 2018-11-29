from BITS.utils import run_command, sge_nize


def run_mapping(contig_fasta="contigs.fasta",
                read_fasta="reads.fasta",
                contig_db_prefix="CONTIGS",
                read_db_prefix="READS",
                script_fname="run_damapper.sh",
                job_submit="qsub",
                n_core="1"):
    """
    Execute DAMAPPER and dump the results.
    """

    n_core = int(n_core)

    script = sge_nize('\n'.join([f"fasta2DAM {contig_db_prefix} {contig_fasta}",
                                 f"DBsplit -s500 {contig_db_prefix}",
                                 f"fasta2DB {read_db_prefix} {read_fasta}",
                                 f"DBsplit -s500 {read_db_prefix}",
                                 f"HPC.damapper -T{n_core} -C -N {contig_db_prefix} {read_db_prefix} | bash -v",
                                 f"LAcat {contig_db_prefix}.{read_db_prefix}.*.las > {contig_db_prefix}.{read_db_prefix}.las",
                                 f"LAsort -a {contig_db_prefix}.{read_db_prefix}.las",
                                 f"LAdump -c {contig_db_prefix} {read_db_prefix} {contig_db_prefix}.{read_db_prefix}.S.las > {contig_db_prefix}.{read_db_prefix}.ladump",
                                 f"DBdump -r -h {contig_db_prefix} > {contig_db_prefix}.dbdump",
                                 f"DBdump -r -h {read_db_prefix} > {read_db_prefix}.dbdump"]),
                      job_name="run_damapper",
                      n_core=n_core,
                      wait=False)

    with open(script_fname, 'w') as f:
        f.write(script)

    run_command(f"{job_submit} {script_fname}")   # wait until finish

