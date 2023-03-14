from os import remove
from os.path import abspath, dirname, join

import pandas as pd
from Bio import AlignIO, SeqIO
from Bio.Align import MultipleSeqAlignment as MSA
from Bio.Seq import Seq


class Slice:
    def __init__(self, input_file, folder="workdir"):
        """Create a Slice object with file paths"""

        # define the file paths
        dir = abspath(join(dirname(__file__), "..", folder))
        self.input_file = input_file
        self.input_path = join(dir, input_file)
        self.output_path = join(dir, input_file[:-6] + "_sliced.fasta")

        # remove the file if it already exists
        try:
            remove(self.output_path)
        except:
            pass

    def slice(self):
        """Function to slice aligned sequences if all sequences have the same aa in the same position"""

        # read alignment file
        seq_records = list(SeqIO.parse(self.input_path, "fasta"))
        len_before = len(seq_records[0].seq)

        # create a dataframe with the sequences to slice easier
        df = pd.DataFrame([list(record.seq) for record in seq_records])

        # get columns to drop
        cols_to_drop = df.nunique()[df.nunique() == 1].index

        # drop the columns
        df.drop(cols_to_drop, axis=1, inplace=True)

        # change sequence in seq record to sliced seq
        for index, seq in df.iterrows():
            seq_sliced = Seq("".join(seq))
            seq_records[index].seq = seq_sliced

        # write msa (multiple sequences alignment) in a file
        AlignIO.write(MSA(seq_records), self.output_path, "fasta")

        len_after = len(seq_records[0].seq)
        ratio = "{:.2f}".format(((len_before - len_after) / len_before) * 100)
        print(f"\nSLICED = Sequences size: {len_before} -> {len_after} (-{ratio}%)")
