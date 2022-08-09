# -*- coding: utf-8 -*-

"""
1. Mutate FVIII human amino acids sequence with a point mutation.
2. Mutate several amino acid sequences separately and align them.
3. Generate alignment file with all mutated sequences aligned.
4. Generate sliced alignment file eliminating repeated amino acids in the same position for all sequences.
"""
from auxiliar import exception_handler, print_finished_time, read_pm_file, get_name_description_seq_file, colorize_svg

from Bio.Data.IUPACData import protein_letters_3to1 as prot3to1
from Bio.Align import MultipleSeqAlignment as MSA
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq, MutableSeq
from Bio import SeqIO, AlignIO

from os.path import abspath, join, dirname
from numpy import NaN
import pandas as pd

import time


# --------------- THIS IS IMMUTABLE during execution --------------------------- #


# diretories
input_dir = abspath(join(dirname(__file__), '..', 'datasets'))
output_dir = abspath(join(dirname(__file__), '..', 'workdir'))

# input files
wild_input_file = 'Human_FVIII_prot.fasta'
wild_input_path = join(input_dir, wild_input_file)

pm_input_file = 'FVIII_point_mutations_v1.csv'
pm_input_path = join(input_dir, pm_input_file)

name, description = get_name_description_seq_file(wild_input_path)

# dataframe with point mutation database cleaned
df_original = read_pm_file(pm_input_path)
df = df_original.copy()
print(df_original)

# ----------------------------- SLICE -------------------------------------------- #


@exception_handler
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

# ----------------------------- ALIGN -------------------------------------------- #


@exception_handler
def align():
    '''function to align all mutated sequences and write in a file'''
    # groupby effect to align sequences with same effect
    df_groupby_effect = df.groupby('effect')

    list_align_path = []
    for effect, df_effect in df_groupby_effect:
        # naming path to write alignment file
        align_file = f'{wild_input_file[:16]}_{effect.lower()}_mutated_{len(df_effect)}_seqs_aligned.fasta'
        align_path = join(output_dir, align_file)

        # list of sequence records to align
        def list_seq_records():
            '''function to create a list of sequence records with all mutated sequences'''
            seq_records = []

            def seq_record(row):
                '''function to create a sequence record using mutated sequence, mutation id, name and description'''
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
        AlignIO.write(mutated_msa, align_path, 'fasta')

        # append align file path to slice later
        list_align_path.append(align_path)

    # return only the missense alignment path
    for align_path in list_align_path:
        if 'missense' in align_path:
            return mutated_msa, align_path

# ----------------------------- MUTATE ------------------------------------------- #


@ exception_handler
def mutate():
    '''function to point mutate all sequences in the dataframe'''
    def mutate_seq(row):
        '''function to point mutate a wild sequence'''
        # get wild sequence and make it mutable
        wild_sequence = MutableSeq(SeqIO.read(wild_input_path, 'fasta').seq)

        # index of python is the position of the mutation - 1
        index = row.position_hgvs - 1

        # amino acids capitalized
        wild = row.wild_aa
        new = row.new_aa

        # verify if the aminoacid i want to mutate in the position is the same in wild sequence
        can_mutate = prot3to1[wild] == wild_sequence[index]

        if can_mutate:
            effect = row.effect.strip().lower()
            # missense just change one amino acid in a specific position
            if effect == 'missense':
                wild_sequence[index] = prot3to1[new]
            # nonsense add a stop codon that eliminate the rest of the sequence
            elif effect == 'nonsense':
                wild_sequence[index:] = ''

            return str(wild_sequence)

        else:
            return NaN

    def mutation_id(row):
        '''function to create an id for the mutation'''
        new = row.new_aa
        wild = row.wild_aa
        sev = row.severity[:3]
        hgvs = row.position_hgvs
        id = f'{sev}_p{{{wild}{hgvs}{new}}}'

        return id

    # create id for each mutated sequence
    df['mutated_seq_id'] = df.apply(mutation_id, axis=1)

    # mutate all sequences
    df['mutated_seq'] = df.apply(mutate_seq, axis=1)

    # drop rows with NaN value (could not mutate)
    print(f'Cannot mutate {df.isna().sum().sum()} sequences.')
    df.dropna(inplace=True)

    # pad sequences using the max sequence length
    max_len = max(len(s) for s in df.mutated_seq.values)
    df['mutated_seq_padded'] = df.mutated_seq.apply(
        lambda s: s.ljust(max_len, 'x'))

    return df


# ----------------------------- MAIN --------------------------------------------- #

def sample_df(percentage=1):
    '''get a sample of the dataframe based on the percentage'''

    # get a sample from the dataframe
    df_sample = df.sample(frac=percentage, random_state=1).sort_index()

    # text to detail database
    database_text = f'\n-- Using {int(percentage*100)}% of the database: {pm_input_path}.'
    mutation_text = f'-- There are {len(df_sample)} point mutations - same number of possible mutated sequences.\n'
    print(database_text)
    print(mutation_text)

    return sample_df


def compare_slicing(seq_before, seq_after):
    '''function to compare msa before and after slice'''

    len_before = len(seq_before[0].seq)
    len_after = len(seq_after[0].seq)
    ratio = '{:.2f}'.format(((len_before - len_after) / len_before)*100)

    print(f'\nSequences length comparison\n')
    print(f'- BEFORE Slice: {len_before} amino acids\n')
    print(f'- AFTER Slice: {len_after} amino acids\n')
    print(f'- Reduced sequence length in {ratio}%\n')


if __name__ == '__main__':
    start_time = time.time()

    # 1. Mutate sequences
    mutated_df = mutate()
    # 2. Align mutated sequences
    aligned_seq_records, missense_align_path = align()
    # 3. Slice mutated sequences to decrease size
    sliced_seq_records = slice(missense_align_path)
    # 4. Sliced sequences length comparison
    compare_slicing(aligned_seq_records, sliced_seq_records)

    print_finished_time(start_time, 'Mutation | Alignment | Slice')
