# -*- coding: utf-8 -*-
import time

from colorize_svg import ColorizeSVG
from point_mutation import PointMutation
from slice import slice, compare_slicing


if __name__ == '__main__':
    start_time = time.time()

    # 0. Testing colorization of phylogenetic tree from seaview
    colorize = ColorizeSVG('results')
    colorize.colorize_seaview_phylo_tree()

    # 1. Read and Clean 'FVIII_point_mutations_v1.csv'
    pm_fviii = PointMutation()
    print(pm_fviii.df.head())

    # 2. Mutate the wild Human FVIII protein using point mutation
    pm_fviii.mutate_sequences()
    print(pm_fviii.df.head())

    # 3. Align all mutated Human FVIII and generate a fasta file
    aligned_seq_records, missense_align_path = pm_fviii.align_mutated_sequences()
    print(pm_fviii.df.head())

    # 4. Use file created on step 3 to slice all mutated Human FVIII and gererate a fasta file
    sliced_seq_records = slice(missense_align_path)
    compare_slicing(aligned_seq_records, sliced_seq_records)
