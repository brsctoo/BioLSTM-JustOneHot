"""
Docstring for validation
"""

import os

from genbank_searcher import BASE_DIR
import modeling
import numpy as np
import pickle
import random
import keras
from minineedle import needle, smith, core

def smooth_predict(predict_raw, window_size=6):
    half = window_size // 2
    smoothed = []

    for i in range(len(predict_raw)):
        start = max(0, i - half)
        end = min(len(predict_raw), i + half + 1)

        window = predict_raw[start:end]

        if sum(window) >= len(window) / 2:
            smoothed.append(1)
        else:
            smoothed.append(0)

    return smoothed

def validate_model(model_path, data_test):
    model = keras.models.load_model(model_path)
    metricas_finais = []

    data_test = pickle.load(open(data_test, "rb"))

    count = 0
    for sample in data_test:
        count +=1
        print("Sequencia:", count, "Progresso:", 100*count/len(data_test), "%")

        tagged_sequence = modeling.tag_positions(sample) # returns something like: [0,0,1,0,0...]
        windows = modeling.slide_window(sample) # sliding window -> list of windows

        X = []
        Y = []

        for j in range(len(windows)):
            X.append(windows[j])
            Y.append(tagged_sequence[j])

        X = np.array(X)
        Y = np.array(Y)

        predict_raw = (model.predict(X) > 0.5).astype("int32") # type: ignore
        predict_suavized = smooth_predict(predict_raw, window_size=6)

        # Final sequence of predicted exons and introns
        final_seq = []
        for i in range(len(predict_suavized)):
            if predict_suavized[i] == 1:
                final_seq.append(sample["sequence"][i])

        # Sequence of true exons and introns
        true_final_seq = []
        for i in range(len(Y)):
            if Y[i] == 1:
                true_final_seq.append(sample["sequence"][i])

        if(len(final_seq) != 0 and len(true_final_seq) != 0 ):
            aligment = needle.NeedlemanWunsch("".join(final_seq), "".join(true_final_seq)) # type: ignore
            aligment.align()

            # 3. Pegamos as sequências alinhadas (com os gaps inseridos)
            seq1_aligned, seq2_aligned = aligment.get_aligned_sequences()

            # 4. Calculamos a Identidade (Matches / Tamanho Total do Alinhamento)
            matches = sum(1 for a, b in zip(seq1_aligned, seq2_aligned) if a == b)
            identidade = matches / len(seq1_aligned) # Retorna um valor perfeito entre 0.0 e 1.0

            print("Identidade de Alinhamento:", identidade)
            metricas_finais.append([identidade])

    print(metricas_finais)

    q0,q1,q2,q3,q4,q5,q6,q7,q8,q9 = 0,0,0,0,0,0,0,0,0,0
    total = 0
    for metrica in metricas_finais:
        # Removemos o abs(). O valor agora é a Identidade real (ex: 0.85 para 85% de match)
        valor_identidade = metrica[0]
        total += valor_identidade

        # Agora essa separação fará sentido real para avaliar o modelo!
        if valor_identidade <= 0.1:
            q0 += 1
        elif valor_identidade <= 0.2:
            q1 += 1
        elif valor_identidade <= 0.3:
            q2 += 1
        elif valor_identidade <= 0.4:
            q3 += 1
        elif valor_identidade <= 0.5:
            q4 += 1
        elif valor_identidade <= 0.6:
            q5 += 1
        elif valor_identidade <= 0.7:
            q6 += 1
        elif valor_identidade <= 0.8:
            q7 += 1
        elif valor_identidade <= 0.9:
            q8 += 1
        else:
            q9 += 1

    media = total / len(metricas_finais)

    print("\n--- RESULTADOS GERAIS ---")
    print("Total de amostras: ", len(metricas_finais))
    print("IDENTIDADE MÉDIA DO ALINHAMENTO: {:.2f}%".format(media * 100))
    print(" 0-10%:   ", q0)
    print(" 10-20%:  ", q1)
    print(" 20-30%:  ", q2)
    print(" 30-40%:  ", q3)
    print(" 40-50%:  ", q4)
    print(" 50-60%:  ", q5)
    print(" 60-70%:  ", q6)
    print(" 70-80%:  ", q7)
    print(" 80-90%:  ", q8)
    print(" 90-100%: ", q9)
    print("\n")
