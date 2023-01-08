# -*- coding: utf-8 -*-
from Bio.Align import MultipleSeqAlignment as MSA
from Bio.Seq import Seq
from Bio import SeqIO, AlignIO

import pandas as pd


def slice(align_path):
    '''function to slice aligned sequences if all sequences have the same aa in the same position'''

    # defining slice output path
    sliced_path = align_path[:-6] + '_sliced.fasta'

    # read alignment file
    aligned_seq_records = list(SeqIO.parse(align_path, 'fasta'))

    # create a dataframe with the sequences to slice easier
    list_sequences = [list(record.seq) for record in aligned_seq_records]
    df_sequences = pd.DataFrame(list_sequences)

    # drop columns that have same amino acid
    nunique = df_sequences.nunique()
    cols_to_drop = nunique[nunique == 1].index
    df_sequences = df_sequences.drop(cols_to_drop, axis=1)

    # change seq in seq record to sliced seq
    for index, prot in df_sequences.iterrows():
        sliced_sequence = Seq(''.join(prot))
        aligned_seq_records[index].seq = sliced_sequence

    # msa = multiple sequences alignment
    sliced_msa = MSA(aligned_seq_records)

    # write msa in a file
    AlignIO.write(sliced_msa, sliced_path, 'fasta')

    return aligned_seq_records


def compare_slicing(seq_before, seq_after):
    '''function to compare msa before and after slice'''

    len_before = len(seq_before[0].seq)
    len_after = len(seq_after[0].seq)
    ratio = '{:.2f}'.format(((len_before - len_after) / len_before)*100)

    print(f'\nSequences length comparison\n')
    print(f'- BEFORE Slice: {len_before} amino acids\n')
    print(f'- AFTER Slice: {len_after} amino acids\n')
    print(f'- Reduced sequence length in {ratio}%\n')
