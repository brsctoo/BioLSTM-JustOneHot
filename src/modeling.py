"""
The script processes genomic sequences with annotated exons and introns to prepare training data.

Steps involved:
1. Load data: Reads a pickled list of sequences and their exon coordinates.
2. Tag positions: Converts each sequence into a list of labels (0 = exon, 1 = intron) of the same length as the sequence.
3. Encode nucleotides: Transforms each base into a numeric representation (A=1, T=2, G=3, C=4, degenerate=5).
4. Create sliding windows: For every position in each sequence, generates a centered window of 60 nucleotides.
5. Pair input with label: Stores each window along with the central position label (0 or 1) in a list XY.
6. Count exon/intron positions: Computes total numbers of introns and exons.
7. Save processed data: Pickles the final XY dataset for later use in model training.

obs.
RYSWKMBDHVN = degenerate bases

sample = {
            "sequence": seq [ACTG...],
            "exon_intervals": exons_intervals [(start1, end1), (start2, end2), ...],
            "exons": exons [ACTG..., ACTG..., ...],
            "intron_intervals": introns_intervals [(start1, end1), (start2, end2), ...],
            "introns": introns [ACTG..., ACTG..., ...],
}

data.append(sample)
"""

import numpy as np
import os
import pickle

# Dictionary to show the one hot enconding arrays
# The normalized values are the probability of the base (key) to be wich nucleotide
# Ex.:
# R gave 50% chance of beeing A and 50% chance of beeing G -> [0.5, 0.0, 0.5, 0.0]
#
# obs:
# N is a unknown base, this means that it has the same chance of beeing any of the them -> [0.25, 0.25, 0.25, 0.25]

BASE_TO_VECTOR = {
    'A':[1,0,0,0],
    'T':[0,1,0,0],
    'G':[0,0,1,0],
    'C':[0,0,0,1],

    'R':[0.5,0,0.5,0],
    'Y':[0,0.5,0,0.5],
    'S':[0,0,0.5,0.5],
    'W':[0.5,0.5,0,0],
    'K':[0,0.5,0.5,0],
    'M':[0.5,0,0,0.5],

    'B':[0,1/3,1/3,1/3],
    'D':[1/3,1/3,1/3,0],
    'H':[1/3,1/3,0,1/3],
    'V':[1/3,0,1/3,1/3],

    'N':[0.25,0.25,0.25,0.25]
}

# vars to cout the qqu
quantidade_bases_degeneradas = 0
quantidade_bases = 0

def tag_positions(sample) -> list[int]:
    """Tag each position in the sequence as exon (1) or intron (0)."""

    tag = [0] * len(sample["sequence"])  # Initialize all positions as intron (0)

    for start, end in sample["exon_intervals"]:
        for i in range(start, end + 1): # Inclusive end position
            tag[i] = 1  # Mark exon positions as 1

    return tag

def slide_window(sample, window_size=60) -> list[list[int]]:
    """Create a sliding window centered at the given position."""

    half = window_size // 2
    seq = transform_baseSeq_to_onehot(sample["sequence"])

    windows = []

    for k in range(len(seq)):
        window = []

        for offset in range(-half, half):
            pos = k + offset
            if pos < 0 or pos >= len(seq):
                window.append([0,0,0,0])  # Padding
            else:
                window.append(seq[pos])

        windows.append(window)

    return windows

def transform_baseSeq_to_onehot(baseSeq) -> list[int]:
    """Transform nucleotide base sequence to numeric sequence."""

    encoded = []
    global quantidade_bases
    quantidade_bases += len(baseSeq)

    for base in baseSeq.upper():
        if base in BASE_TO_VECTOR:
            if base not in ['A', 'T', 'G', 'C']:
                global quantidade_bases_degeneradas
                quantidade_bases_degeneradas += 1

        encoded.append(
            BASE_TO_VECTOR.get(base, [0,0,0,0])  # unknown → padding
        )

    return np.array(encoded, dtype=np.float16) # type: ignore

def save_XY_to_file(output_path, X_list, y_list):
    """
    Recebe listas separadas para evitar o uso de np.array(dtype=object).
    """
    print(f"Empilhando {len(X_list)} amostras... Isso pode demorar um pouco.")

    # Converte a lista de janelas diretamente para um array 4D
    # X shape original: (N, 60, 4)
    X = np.stack(X_list).astype(np.float16)

    # Converte labels para int8 (ocupa muito menos espaço que o padrão int64)
    y = np.array(y_list).astype(np.int8)

    print(f"Salvando arquivo comprimido em: {output_path}")
    np.savez_compressed(output_path, X=X, y=y)

    # Limpeza explícita para garantir liberação de RAM
    del X
    del y

def build_XY_dataset(data):
    """Build the X and y lists separately to avoid object arrays."""
    X_final = []
    y_final = []

    for sample in data:
        print("Processing: ", sample["sequence"][:50], "...")
        tagged_seq = tag_positions(sample)
        windows = slide_window(sample) # windows já retorna np.float16

        # Em vez de XY.append([window, label]), separamos em duas listas
        for j in range(len(windows)):
            X_final.append(windows[j])
            y_final.append(tagged_seq[j])

    return X_final, y_final

def modeling_train_data(data_filepath_input, XY_filepath_output):
    data = pickle.load(open(data_filepath_input, "rb"))

    # Captura as listas processadas
    X_list, y_list = build_XY_dataset(data)

    # Save using
    save_XY_to_file(XY_filepath_output, X_list, y_list)

    # Debug of the values
    print("Concluído com sucesso!")
    print("Quantidade de bases degeneradas: ", quantidade_bases_degeneradas)
    print("Quantidade de bases: ", quantidade_bases)
    ratio_degenerate = (quantidade_bases_degeneradas / quantidade_bases) * 100
    print("Porcentagem de bases degeneradas: ", ratio_degenerate, "%")

    return ratio_degenerate
