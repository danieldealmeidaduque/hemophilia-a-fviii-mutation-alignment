from mutate import PointMutateFVIII

if __name__ == "__main__":

    mutate = PointMutateFVIII(
        input_wild="Human_FVIII_prot.fasta",
        input_champ="champ-mutation-list-q4-clean.xlsx",
        input_folder="datasets",
        output_folder="workdir",
    )

    mutate.mutate_sequences()
