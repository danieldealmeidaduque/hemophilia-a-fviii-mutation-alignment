from os.path import abspath, join, dirname
from os import remove
from Bio import SeqIO

import pandas as pd
import time
import sys
import re

# ------------------------ Auxiliar Functions ------------------------------------ #


def exception_handler(func):
    '''cleaner way to handle exceptions and exit the program'''
    def inner_function(*args, **kwargs):
        try:
            retry = func(*args, **kwargs)
            return retry
        except:
            print(f'\n\t !!!!!!! ERROR IN FUNCTION: {func.__name__}\n')
            sys.exit(1)
    return inner_function


@ exception_handler
def print_finished_time(start_time, finish_msg):
    '''print finished time formatted'''
    f_time = time.time()
    e_time = round((f_time - start_time), 2)
    msg = f'\n\t{finish_msg} - Finished  in {e_time} s.\n'
    print(msg)


@ exception_handler
def print_highlighted(s):
    '''print upper and lower highlighted text'''
    print('\n\t\t\t ----------- ' + s.upper() + ' ------------\n')


def strip_string(s):
    try:
        s = s.strip()
    finally:
        return s


# ------------------------ Read Input Files -------------------------------------- #


@exception_handler
def get_name_description_seq_file(input_path):
    '''read the wild protein fasta file to get the sequence information'''
    seq_record = SeqIO.read(input_path, 'fasta')

    name = seq_record.name
    description = seq_record.description

    return name, description


@ exception_handler
def read_pm_file(input_path, csv_sep='\t'):
    '''read and clean the point mutations input file to a dataframe'''
    df = pd.read_csv(input_path, sep=csv_sep)
    # remove extra spaces from string values
    df = df.applymap(lambda s: strip_string(s))
    # padronize column names to access with dot
    df.columns = [c.lower().replace(' ', '_') for c in df.columns]
    # set case_id as index
    df = df.set_index('case_id')

    # func to clean protein_change column to clear view
    def clean_prot_change(prot_change):
        pattern = '(?<=\()(.*?)(?=\))'
        clean_str = re.search(pattern, prot_change)
        clean_str = clean_str.group().replace(' ', '')
        clean_str = clean_str.strip()  # .lower()
        return clean_str

    # func to get only the wild amino acid from protein_change column
    def get_wild_aa(prot_change):
        return prot_change[:3]

    # func to get only the new amino acid from protein_change column
    def get_new_aa(row):
        effect = row.effect.strip()
        if effect == 'Missense':
            return row.protein_change[-3:]
        else:
            return '---'

    # apply changes in the dataframe
    df.protein_change = df.protein_change.apply(clean_prot_change)
    df['wild_aa'] = df.protein_change.apply(get_wild_aa)
    df['new_aa'] = df.apply(get_new_aa, axis=1)

    return df


@ exception_handler
def colorize_svg(input_path):
    # diretory
    dir = abspath(join(dirname(__file__), '..', 'workdir'))

    # input
    input_file = 'phyml_Blosum62_test1.svg'
    input_path = join(dir, input_file)

    # output
    output_file = input_file[:-4] + '_colorized.svg'
    output_path = join(dir, output_file)

    # if output file already exists then delete it
    try:
        remove(output_path)
    except:
        print('SVG colorido NÃO existe. Criando.....')

    '''Colorize the id's of the phylogenetic tree SVG file from the seaview'''
    with open(input_path) as f1, open(output_path, 'a') as f2:
        for line in f1:
            if 'Sev' in line:
                color = 'red'
            elif 'Mod' in line:
                color = 'blue'
            elif 'Mil' in line:
                color = 'green'
            else:
                color = 'black'

            line = line.replace('rgb(0,0,0)', color)
            f2.write(line)
