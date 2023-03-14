from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment as MSA
from os.path import abspath, dirname, join
from mutate import Mutate


class Align:
    """Class representing the point mutations of the csv file as a dataframe"""

    def __init__(
        self, input_file, output_folder="workdir", input_folder="datasets", sep="\t"
    ):
        """Create a Align object"""

        # define the file paths
        self.input_dir = abspath(join(dirname(__file__), "..", input_folder))
        self.output_dir = abspath(join(dirname(__file__), "..", output_folder))
        self.input_file = input_file
        self.input_path = join(self.input_dir, input_file)

        # create dataframe using csv file
        self.df = Mutate("Human_FVIII_prot.fasta", input_file)

    def align_mutated_sequences(self):

        seq_record = SeqIO.read(self.input_path, "fasta")
        name = seq_record.name
        description = seq_record.description

        # groupby effect to align sequences with same effect
        df_groupby_effect = self.df.groupby("effect")

        list_align_path = []
        for effect, df_effect in df_groupby_effect:
            # naming path to write alignment file
            align_file = f"{self.input_file[:16]}_{effect.lower()}_mutated_{len(df_effect)}_seqs_aligned.fasta"
            align_path = join(self.output_dir, align_file)

            # list of sequence records to align
            def list_seq_records():
                """function to create a list of sequence records with all mutated sequences"""
                seq_records = []

                def seq_record(row):
                    """function to create a sequence record using mutated sequence, mutation id, name and description"""
                    seq = Seq(row.mutated_seq_padded)
                    id = row.mutated_seq_id
                    seqRecord = SeqRecord(seq, id, name, description)

                    seq_records.append(seqRecord)

                df_effect.apply(seq_record, axis=1)

                return seq_records

            # get list of all mutated sequence records using the dataframe
            mutated_seq_records_to_align = list_seq_records()

            # msa = multiple sequences alignment
            mutated_msa = MSA(mutated_seq_records_to_align)

            # write msa in a file
            AlignIO.write(mutated_msa, align_path, "fasta")

            # append align file path to slice later
            list_align_path.append(align_path)

        # return only the missense alignment path
        for align_path in list_align_path:
            if "missense" in align_path:
                return mutated_msa, align_path
