from mutate import Mutate

if __name__ == "__main__":

    mutate = Mutate(
        input_wild="Human_FVIII_prot.fasta",
        input_champ="champ-mutation-list-q4-clean.xlsx",
        folder="datasets",
    )