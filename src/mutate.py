import re
from os.path import abspath, dirname, join

import pandas as pd
from Bio import AlignIO, SeqIO
from Bio.Align import MultipleSeqAlignment as MSA
from Bio.Data.IUPACData import protein_letters_3to1 as prot3to1
from Bio.Seq import MutableSeq, Seq
from Bio.SeqRecord import SeqRecord
from numpy import NaN


def get_info(s):
    wild = re.search("^\D+", s).group()
    pos = re.search("\d+", s).group()
    new = re.search("\D+$", s).group()
    return wild, int(pos), new


class Mutate:
    """Class representing the point mutations of the csv file as a dataframe"""

    def __init__(self, input_wild, input_champ, folder):
        """Create a Mutate object"""

        # define the file paths
        dir = abspath(join(dirname(__file__), "..", folder))
        self.path_champ = join(dir, input_champ)
        self.path_wild = join(dir, input_wild)

        # create dataframe using the champ mutation list in excel
        df = pd.read_excel(self.path_champ)

        # create a dict with severity as key and list of mutations as values
        self.d = {}
        for sev, prot in df.groupby("Reported Severity"):
            self.d[sev] = prot["HGVS Protein"].to_list()

        # mutate all sequences in the dataframe
        self.mutate_sequences()

        # self.df.to_excel("teste.xlsx")

    def mutate_sequences(self):
        """Apply all mutations from the champ mutation list"""

        def mutate_seq(prot):
            """Mutate a single sequence"""

            # get the wild amino acid sequence and make it mutable
            wild_seq = MutableSeq(SeqIO.read(self.path_wild, "fasta").seq)

            # alias to faciliate manipulation
            wild, pos, new = get_info(prot)

            # verify if wild is a valid amino acid
            if not (wild in prot3to1.keys()):
                return NaN

            # verify if wild is the same amino acid in the wild sequence to mutate it
            if not (prot3to1[wild] == wild_seq[pos - 1]):
                return NaN

            # * is stop codon
            if new == "*":
                wild_seq[pos - 1 :] = ""
            else:
                wild_seq[pos - 1] = prot3to1[new]

            return str(wild_seq)

        seq_record = SeqIO.read(self.path_wild, "fasta")
        name = seq_record.name
        description = seq_record.description

        seq_records = []
        for sev, prots in self.d.items():
            for prot in prots:
                id = sev + "{" + prot + "}"

                seq = mutate_seq(prot)
                if seq is NaN:
                    continue
                else:
                    seq_record = SeqRecord(MutableSeq(seq), id, name, description)
                    seq_records.append(seq_record)

        # pad all mutated sequences
        max_len = max(len(str(r.seq)) for r in seq_records)
        print(max_len)

        for r in seq_records:
            s = str(r.seq).ljust(max_len, "x")
            r.seq = Seq(s)

        # generate fasta file with mutated sequence records
        AlignIO.write(MSA(seq_records), f"align_path{len(seq_records)}.fasta", "fasta")


mutate = Mutate(
    input_wild="Human_FVIII_prot.fasta",
    input_champ="champ-mutation-list-q4-clean.xlsx",
    folder="datasets",
)
