import pandas as pd
from Bio import SeqIO
from Bio.Data.IUPACData import protein_letters_3to1 as prot3to1
from Bio.Seq import MutableSeq
from os.path import abspath, dirname, join
from numpy import NaN


class Mutate:
    """Class representing the point mutations of the csv file as a dataframe"""

    def __init__(self, input_wild, input_file, folder="datasets", sep="\t"):
        """Create a Mutate object"""

        # define the file paths
        input_dir = abspath(join(dirname(__file__), "..", folder))
        self.input_file = input_file
        self.input_path = join(input_dir, input_file)
        self.input_wild_path = join(input_dir, input_wild)

        # create dataframe using csv file
        self.df = pd.read_csv(self.input_path, sep=sep)

        # mutate all sequences in the dataframe
        self.mutate_sequences()

    def mutate_sequences(self):
        """Apply all mutations of the dataframe in the wild sequence"""

        def mutate_seq(row):
            """Mutate a single sequence"""

            # get the wild sequence and make it mutable
            wild_seq = MutableSeq(SeqIO.read(self.input_wild_path, "fasta").seq)

            # get alias for used informations
            index = row.mut_pos - 1
            wild_aa = row.wild_aa.strip()
            new_aa = row.new_aa.strip()

            # verify if the wild amino acid is the same as the original
            can_mutate = prot3to1[wild_aa] == wild_seq[index]

            if not can_mutate:
                return NaN
            else:
                # missense just change one amino acid in a specific position
                if row.effect == "Missense":
                    wild_seq[index] = prot3to1[new_aa]

                # nonsense add a stop codon that eliminate the rest of the sequence
                elif row.effect == "Nonsense":
                    wild_seq[index:] = ""

                return str(wild_seq)

        def mutate_id(row):
            """Create an id for the mutated sequence"""
            new_aa = row.new_aa
            wild_aa = row.wild_aa
            mut_pos = row.mut_pos
            sev = row.severity[:3]
            return f"{sev}_p{{{wild_aa}{mut_pos}{new_aa}}}"

        # alias to improve readability
        df = self.df

        # create an id for each mutated sequence
        df["mutated_id"] = df.apply(mutate_id, axis=1)

        # apply the mutation on each sequence
        df["mutated_seq"] = df.apply(mutate_seq, axis=1)

        # drop rows that could not mutate
        print(f"Could not mutate {df.isna().sum().sum()} sequences.")
        df.dropna(inplace=True)

        # pad all mutated sequences
        max_len = max(len(s) for s in df.mutated_seq.values)
        df["mutated_seq"] = df.mutated_seq.apply(lambda s: s.ljust(max_len, "x"))
