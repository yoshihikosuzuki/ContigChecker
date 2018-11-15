from typing import List
from dataclasses import dataclass, field
import pandas as pd
from BITS.utils import run_command, sge_nize

script_fname = "run_damapper.sh"
job_submit = "qsuball"


def load_dbdump(dbdump_fname):
    """
    Load DAZZ_DB ID, fasta header name, and sequence length for each read in a database.
    """

    seqs = {}
    with open(dbdump_fname, 'r') as f:
        for line in f:
            data = line.strip().split(' ')
            if data[0] == "R":
                dbid = int(data[1])
            elif data[0] == "H":
                header = data[2][1:]
            elif data[0] == "L":   
                length = int(data[3])
                seqs[dbid] = [header, length]

    seqs = pd.DataFrame.from_dict(seqs, orient="index", columns=("header", "length"))
    seqs.index.name = "dbid"
    return seqs


@dataclass(repr=False, eq=False)
class Mapping:
    """
    Class of a mapping of a read onto a contig.
    Mapping type is also determined based on the following features:

    * Number of split alignments: single or chain
    * Terminal                  : proper or clipped
    * Inclusion relationship    : contains, contained, start_dovetail, or end_dovetail

    Note that "contains" (a read contains a contig from end to end) can occur only when terminal is "proper".
    """

    # Initialization = Loading a "P" line in a ladump file
    contig_id: int
    read_id: int
    contig_len: int
    read_len: int
    strand: str
    split_type: str = field(default="single", init=False)   # "single" if not chain else "chain"
    break_starts: List[int] = field(default=[], init=False)   # intervals of breaks between chains
    break_ends: List[int] = field(default=[], init=False)

    def add_alignment(self, chain, contig_s, contig_e, read_s, read_e):
        """
        Add an alignment (which might be part of a chain) = Loading a "C" line in a ladump file.
        """

        if not chain:
            self.contig_s = contig_s
            self.contig_e = contig_e
            self.read_s = read_s
            self.read_e = read_e
        else:
            self.split_type = "chain"
            self.break_starts.append(self.contig_e - 1)   # TODO: 1 is necessary?
            self.break_ends.append(contig_s)
            self.contig_e = contig_e   # Update only end positions
            self.read_e = read_e
            
    def _check_type(self):
        if self.contig_s == 0 and self.contig_e == self.contig_len:
            self.terminal = "proper"
            self.map_type = "contains"
        elif self.contig_s == 0 and self.read_s > 0:
            self.map_type = "start_dovetail"
            self.terminal = "proper" if self.read_e == self.read_len else "clipped"
        elif self.contig_e == self.contig_len and self.read_e < self.read_len:
            self.map_type = "end_dovetail"
            self.terminal = "proper" if self.read_s == 0 else "clipped"
        else:
            self.map_type = "contained"
            self.terminal = "proper" if self.read_s == 0 and self.read_e == self.read_len else "clipped"
        
    def as_list(self):
        self._check_type()
        return [self.contig_id,
                self.read_id,
                self.strand,
                self.contig_s,
                self.contig_e - 1,
                self.split_type,
                self.terminal,
                self.map_type,
                self.break_starts,
                self.break_ends]


def load_ladump(ladump_fname, contigs, reads):
    mappings = {}
    index = 0
    with open(ladump_fname, 'r') as f:
        flag_first = True
        for line in f:
            data = line.strip().split(' ')
            if data[0] == "P":
                chain = True if data[4] == "-" else False
                if chain:
                    continue
                
                if flag_first:
                    flag_first = False
                else:
                    # Fix a mapping after all split alignments in a chain are added
                    mappings[index] = mapping.as_list()
                    index += 1
                
                # Prepare a new mapping
                contig_id, read_id, strand = int(data[1]), int(data[2]), data[3]
                mapping = Mapping(contig_id,
                                  read_id,
                                  contigs.loc[contig_id, "length"],
                                  reads.loc[read_id, "length"],
                                  strand)

            elif data[0] == "C":
                contig_s, contig_e, read_s, read_e = list(map(int, data[1:]))
                mapping.add_alignment(chain, contig_s, contig_e, read_s, read_e)

        mappings[index] = mapping.as_list()   # the last mapping

    mappings = pd.DataFrame.from_dict(mappings, orient="index",
                                      columns=("contig_id",
                                               "read_id",
                                               "strand",
                                               "contig_start",
                                               "contig_end",
                                               "split_type",
                                               "terminal_type",
                                               "map_type",
                                               "break_starts",
                                               "break_ends"))
    return mappings


def run_mapping(contigs_fname, reads_fname, n_core):
    # Execute DAMAPPER and dump the results
    script = sge_nize('\n'.join([f"fasta2DAM CONTIGS {contigs_fname}",
                                 f"DBsplit -s500 CONTIGS",
                                 f"fasta2DB READS {reads_fname}",
                                 f"DBsplit -s500 READS",
                                 f"HPC.damapper -T{n_core} -C -N CONTIGS READS | bash -v",
                                 f"LAcat CONTIGS.READS.*.las > CONTIGS.READS.las",
                                 f"LAsort -a CONTIGS.READS.las",
                                 f"LAdump -c CONTIGS READS CONTIGS.READS.S.las > ladump",
                                 f"DBdump -r -h CONTIGS > CONTIGS.dbdump",
                                 f"DBdump -r -h READS > READS.dbdump"]),
                      job_name="run_damapper",
                      n_core=n_core,
                      wait=True)

    with open(script_fname, 'w') as f:
        f.write(script)

    run_command(f"{job_submit} {script_fname}")   # wait until finish

    # Load the dump data and reformat them
    contigs = load_dbdump("CONTIGS.dbdump")
    reads = load_dbdump("READS.dbdump")
    mappings = load_ladump("ladump", contigs, reads)
    #counts = count_start_end(mappings, contigs)

    contigs.to_csv("contigs.csv", sep='\t')
    reads.to_csv("reads.csv", sep='\t')
    mappings.to_csv("mappings.csv", sep='\t')
    #counts.to_csv("counts.csv", sep='\t')
