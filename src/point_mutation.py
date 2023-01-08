# -*- coding: utf-8 -*-
from Bio.Data.IUPACData import protein_letters_3to1 as prot3to1
from Bio.Align import MultipleSeqAlignment as MSA
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq, MutableSeq
from Bio import SeqIO, AlignIO

from os.path import abspath, join, dirname
from numpy import NaN
import pandas as pd
import re


class PointMutation:
    '''Class representing the point mutations of the csv file as a dataframe'''

    def __init__(self, input_folder='datasets', input_file='FVIII_point_mutations_v1.csv', sep='\t'):
        '''Construct a PointMutations object with input path and df attribute'''
        input_dir = abspath(join(dirname(__file__), '..', input_folder))
        input_path = join(input_dir, input_file)
        self.input_path = input_path
        self.df_original = pd.read_csv(input_path, sep=sep)
        self.df = self.df_original.copy()
        self.preprocess_pm_df()

    def get_df_sample(self, percentage=1.0, random_state=None):
        '''Get a sample of the dataframe of the PointMutations'''
        df_sample = self.df.sample(frac=percentage, random_state=random_state)
        print(f'-- Got {int(percentage*100)}% = {len(df_sample)} mutations --')
        return df_sample

    def standardize_index_columns(self):
        '''Standardize columns of the dataframe and set index as 'case_id' column'''
        df = self.df
        # Standardize columns to be able to access with dot notation
        df.columns = [c.lower().replace(' ', '_') for c in df.columns]
        # Set 'case_id' as index of the dataframe
        df.set_index('case_id', inplace=True)

    def strip_string_values(self):
        '''Strip all string values of the dataframe'''
        def strip_string(s):
            try:
                s = s.strip()
            finally:
                return s
        df = self.df
        # Apply the function above to each value of the entire dataframe
        df = df.applymap(lambda s: strip_string(s))

    def create_prot_change_column(self):
        '''Create 'prot_change' columnn that is 'protein change' column after cleannig'''
        def clean_prot_change(prot_change):
            # Get the value inside parentheses
            inside_str = re.search('(?<=\()(.*?)(?=\))', prot_change).group()
            # Return the value inside parentheses stripped
            return inside_str.strip()
        df = self.df
        # Apply the function above to each value of the protein change column
        df['prot_change'] = df.protein_change.apply(clean_prot_change)

    def create_wild_aa_column(self):
        '''Create a column that is the wild amino acid of the mutation'''
        def get_wild_aa(prot_change):
            # Example: 'Phe1696Ser' -> return 'Phe'
            return prot_change[:3]
        df = self.df
        # Apply the function above to each value of the protein change column
        df['wild_aa'] = df.prot_change.apply(get_wild_aa)

    def create_new_aa_column(self):
        '''Create a column that is the new amino acid of the mutation'''
        def get_new_aa(row):
            if row.effect == 'Missense':
                # Example: 'Phe1696Ser' -> return 'Ser'
                return row.prot_change[-3:]
            if row.effect == 'Nonsense':
                # Example: 'Phe1696*' -> return NaN
                return NaN
        df = self.df
        # Apply the function above to each row of the dataframe
        df['new_aa'] = df.apply(get_new_aa, axis=1)

    def preprocess_pm_df(self):
        '''Apply all preprocessing functions'''
        self.standardize_index_columns()
        self.strip_string_values()
        self.create_prot_change_column()
        self.create_wild_aa_column()
        self.create_new_aa_column()

    def mutate_sequences(self, input_folder='datasets', input_file='Human_FVIII_prot.fasta'):
        '''Apply all mutations of the dataframe in the wild sequence'''
        input_dir = abspath(join(dirname(__file__), '..', input_folder))
        input_path = join(input_dir, input_file)

        def mutate_seq(row):
            # Get the wild sequence and make it mutable
            wild_seq = MutableSeq(SeqIO.read(input_path, 'fasta').seq)

            index = row.position_hgvs - 1
            wild_aa = row.wild_aa
            new_aa = row.new_aa

            # Verify If the wild amino acid is the same as the original
            can_mutate = prot3to1[wild_aa] == wild_seq[index]

            if can_mutate:
                # Missense just change one amino acid in a specific position
                if row.effect == 'Missense':
                    wild_seq[index] = prot3to1[new_aa]
                # Nonsense add a stop codon that eliminate the rest of the sequence
                elif row.effect == 'Nonsense':
                    wild_seq[index:] = ''

                return str(wild_seq)
            else:
                return NaN

        def mutation_id(row):
            '''Create an id for the mutation sequence'''
            new = row.new_aa
            wild = row.wild_aa
            sev = row.severity[:3]
            hgvs = row.position_hgvs
            id = f'{sev}_p{{{wild}{hgvs}{new}}}'

            return id

        df = self.df
        # Create id for each mutated sequence
        df['mutated_seq_id'] = df.apply(mutation_id, axis=1)

        # Mutate each sequence
        df['mutated_seq'] = df.apply(mutate_seq, axis=1)

        # Drop rows with NaN value (could not mutate)
        print(f'Cannot mutate {df.isna().sum().sum()} sequences.')
        df.dropna(inplace=True)

        # Pad sequences using the max sequence length
        max_len = max(len(s) for s in df.mutated_seq.values)
        df['mutated_seq_padded'] = df.mutated_seq.apply(
            lambda s: s.ljust(max_len, 'x'))

    def align_mutated_sequences(self, output_folder='workdir', input_folder='datasets', input_file='Human_FVIII_prot.fasta'):
        '''Align all mutated sequences in the dataframe and write in a file'''
        output_dir = abspath(join(dirname(__file__), '..', output_folder))
        input_dir = abspath(join(dirname(__file__), '..', input_folder))
        input_path = join(input_dir, input_file)

        seq_record = SeqIO.read(input_path, 'fasta')
        name = seq_record.name
        description = seq_record.description

        # groupby effect to align sequences with same effect
        df_groupby_effect = self.df.groupby('effect')

        list_align_path = []
        for effect, df_effect in df_groupby_effect:
            # naming path to write alignment file
            align_file = f'{input_file[:16]}_{effect.lower()}_mutated_{len(df_effect)}_seqs_aligned.fasta'
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
