from mutate import PointMutateFVIII

if __name__ == "__main__":

    mutate = PointMutateFVIII(
        input_wild="Human_FVIII_prot.fasta",
        input_champ="champ-mutation-list-q4-clean.xlsx",
        input_folder="data",
        output_folder="result",
    )

    mutate.mutate_sequences()
