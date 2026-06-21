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
DEFAULT_INJECTION_RATE = 0.0
DEFAULT_NAME = "actin_fungi"

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
genbank_data_filepath_input = os.path.join(BASE_DIR, "../assets/genbank_data/actin_fungi")

def get_output_paths(name):
    mod1 = os.path.join(BASE_DIR, f"../assets/processed_data/mod1/data_{name}")
    mod2 = os.path.join(BASE_DIR, f"../assets/processed_data/mod2/data_XY_{name}.npz")
    result = os.path.join(BASE_DIR, f"../assets/result/model_{name}_onehot.h5")
    return mod1, mod2, result

def train_pipeline(injection_rate, name):
    mod1_filepath_output, mod2_filepath_output, result_filepath_output = get_output_paths(name)

    # Use the genbank_data to create the mod1 data, which is the preprocessed data ready for modeling
    print("Starting the training pipeline...")
    print("Injection rate for degenerate nucleotides: ", injection_rate)
    print("Experiment name: ", name)

    # .gb file → mod1 (train and test)
    genbank_reader.save_preprocessed_genbank_file(genbank_data_filepath_input, mod1_filepath_output, injection_rate)

    # Use the mod1 train data to create the mod2 data, which is the X and y arrays ready for model training
    # .mod1 train → mod2 (X and y)
    modeling.modeling_train_data(mod1_filepath_output + "_train.mod1", mod2_filepath_output)

    # Use the mod2 data to train the model and save it for later use in validation.py
    # mod2 (X and y) → h5 model
    train_model.train_model(mod2_filepath_output, result_filepath_output)
    return None

def validate_pipeline(name):
    # Use the model trained and the test data to validate the model and print the final metrics
    mod1_filepath_output, _, result_filepath_output = get_output_paths(name)

    # h5 model + mod1 test → metric
    validation.validate_model(result_filepath_output, mod1_filepath_output + "_test.mod1")

def main():
    parser = argparse.ArgumentParser(description="Pipeline Bi-LSTM para Identificação de Íntrons/Éxons")

    # Argumento opcional: Modo de operação
    parser.add_argument("mode", nargs='?', choices=["train", "test", "full"],
    help="Escolha o modo: train (treinar), test (validar) ou full (ambos)")

    # Argumento opcional: Taxa de injeção de nucleotídeos degenerados
    parser.add_argument("--injection-rate", type=float, default=DEFAULT_INJECTION_RATE,
            help="Taxa de injeção de nucleotídeos degenerados (ex: 0.5 para 50%%). Padrão: 0.0")

    # Argumento opcional: Nome do experimento
    parser.add_argument("--name", type=str, default=DEFAULT_NAME,
        help="Nome do experimento, usado nos arquivos de saída (ex: actin_fungi_rate50). Padrão: actin_fungi")

    args = parser.parse_args()
    injection_rate = args.injection_rate
    name = args.name

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
      print(f"\nTaxa de injeção atual: {injection_rate} (use --injection-rate para alterar)")
      print(f"Nome do experimento: {name} (use --name para alterar)")

      choice = input("\nDigite o número da opção (1/2/3): ").strip()

      modes_map = {"1": "train", "2": "test", "3": "full"}
      args.mode = modes_map.get(choice)

      if args.mode is None:
          print("\n❌ Opção inválida! Use 1, 2 ou 3.")
          return

    if args.mode == "train":
        train_pipeline(injection_rate, name)
    elif args.mode == "test":
        validate_pipeline(name)
    elif args.mode == "full":
        train_pipeline(injection_rate, name)
        validate_pipeline(name)

if __name__ == "__main__":
    main()
