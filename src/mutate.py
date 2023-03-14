import re
from os.path import abspath, dirname, join

import pandas as pd
from Bio import AlignIO, SeqIO
from Bio.Align import MultipleSeqAlignment as MSA
from Bio.Data.IUPACData import protein_letters_3to1 as prot3to1
from Bio.Seq import MutableSeq, Seq
from Bio.SeqRecord import SeqRecord
from numpy import NaN


class PointMutateFVIII:
    def __init__(self, input_wild, input_champ, input_folder, output_folder):
        """Create a Mutate object"""

        # define the file paths
        self.input_dir = abspath(join(dirname(__file__), "..", input_folder))
        self.output_dir = abspath(join(dirname(__file__), "..", output_folder))
        self.path_champ = join(self.input_dir, input_champ)
        self.path_wild = join(self.input_dir, input_wild)
        self.input_wild = input_wild

        # create a dataframe using the champ mutation list
        df = pd.read_excel(self.path_champ)

        # create a dictionary with severity as keys and list of mutations as values
        self.d = {}
        for sev, prot in df.groupby("Reported Severity"):
            self.d[sev] = prot["HGVS Protein"].to_list()

    def mutate_sequences(self):
        """Apply all mutations from the champ mutation list"""

        def get_mutation(s):
            wild = re.search("^\D+", s).group()
            pos = re.search("\d+", s).group()
            new = re.search("\D+$", s).group()
            return wild, int(pos), new

        def mutate_seq(prot):
            """Mutate a single sequence of amino acids"""

            # get the wild amino acid sequence and make it mutable
            wild_seq = MutableSeq(SeqIO.read(self.path_wild, "fasta").seq)

            # alias to faciliate mutation manipulation
            wild_aa, pos, new_aa = get_mutation(prot)

            # verify if wild_aa is a valid amino acid
            if not (wild_aa in prot3to1.keys()):
                return NaN

            # verify if wild_aa is the same amino acid in the wild sequence to be able to mutate
            if not (prot3to1[wild_aa] == wild_seq[pos - 1]):
                return NaN

            # if is stop codon (*) remove all amino acids after it
            if new_aa == "*":
                wild_seq[pos - 1 :] = ""
            else:
                wild_seq[pos - 1] = prot3to1[new_aa]

            # return as string to be able to
            return MutableSeq(wild_seq)

        # get informations from the wild FVIII fasta file
        wild_seq_record = SeqIO.read(self.path_wild, "fasta")
        name = wild_seq_record.name
        description = wild_seq_record.description

        # apply each mutation in a wild FVIII and generate a list of seq records
        seq_records = []
        for sev, prots in self.d.items():
            for prot in prots:
                id = sev + "{" + prot + "}"
                seq = mutate_seq(prot)

                if seq is NaN:
                    continue
                else:
                    seq_record = SeqRecord(seq, id, name, description)
                    seq_records.append(seq_record)

        # pad all mutated sequences to be able to write it into the fasta file
        max_len = max(len(str(r.seq)) for r in seq_records)
        for r in seq_records:
            s = str(r.seq).ljust(max_len, "x")
            r.seq = Seq(s)

        # define output path
        output_file = f"{self.input_wild[:-6]}_{len(seq_records)}_sequences_align.fasta"
        output_path = join(self.output_dir, output_file)

        # generate a fasta file with all mutated sequence records
        AlignIO.write(MSA(seq_records), output_path, "fasta")
