"""
Script that with only one function, pipeline(), that executes the entire pipeline of data loading, preprocessing and model training. 
01 - genbank_reader.py: contains functions to read and preprocess the GenBank file, including validation of records, separation of train and test datasets, and saving the processed data to files.
"""

import os
import numpy as np
import genbank_reader
import modeling
import train_model
import validation

import argparse
# For injection rate, 100% is like 33% of the bases being replaced by degenerate nucleotides, 50% is like 16.5%, and so on.
INJECTION_RATE = 0.01

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
genbank_data_filepath_input = os.path.join(BASE_DIR, "../assets/genbank_data/actin_fungi")
mod1_filepath_output = os.path.join(BASE_DIR, "../assets/processed_data/mod1/data_actin_fungi")
mod2_filepath_output = os.path.join(BASE_DIR,"../assets/processed_data/mod2/data_XY_actin_fungi.npz")
result_filepath_output = os.path.join(BASE_DIR,"../assets/result/model_actin_fungi_onehot.h5")
 
def train_pipeline():
    # Use the genbank_data to create the mod1 data, which is the preprocessed data ready for modeling
    # .gb file → mod1 (train and test)
    print("Starting the training pipeline...")
    print("Injection rate for degenerate nucleotides: ", INJECTION_RATE)
    genbank_reader.save_preprocessed_genbank_file(genbank_data_filepath_input, mod1_filepath_output, INJECTION_RATE)

    # Use the mod1 train data to create the mod2 data, which is the X and y arrays ready for model training
    # .mod1 train → mod2 (X and y)
    modeling.modeling_train_data(mod1_filepath_output + "_train.mod1", mod2_filepath_output)

    # Use the mod2 data to train the model and save it for later use in validation.py
    # mod2 (X and y) → h5 model
    train_model.train_model(mod2_filepath_output, result_filepath_output)

    return None

def validate_pipeline():
    # Use the model trained and the test data to validate the model and print the final metrics
    # h5 model + mod1 test → metric
    validation.validate_model(result_filepath_output, mod1_filepath_output + "_test.mod1")

def main():
    parser = argparse.ArgumentParser(description="Pipeline Bi-LSTM para Identificação de Íntrons/Éxons")

    # Argumento opcional: Modo de operação
    parser.add_argument("mode", nargs='?', choices=["train", "test", "full"], 
    help="Escolha o modo: train (treinar), test (validar) ou full (ambos)")

    args = parser.parse_args()
    
    # Se nenhum modo foi fornecido, exibir menu interativo
    if args.mode is None:
        print("\n" + "="*50)
        print("Pipeline Bi-LSTM para Identificação de Íntrons/Éxons")
        print("="*50)
        print("\nEscolha uma opção:")
        print("1 - train  (Treinar o modelo)")
        print("2 - test   (Validar o modelo)")
        print("3 - full   (Treinar e validar)")
        print("="*50)
        
        choice = input("\nDigite o número da opção (1/2/3): ").strip()
        
        modes_map = {"1": "train", "2": "test", "3": "full"}
        args.mode = modes_map.get(choice)
        
        if args.mode is None:
            print("\n❌ Opção inválida! Use 1, 2 ou 3.")
            return

    if args.mode == "train":
        train_pipeline()
    elif args.mode == "test":
        validate_pipeline()
    elif args.mode == "full":
        train_pipeline()
        validate_pipeline()

if __name__ == "__main__":
    main()