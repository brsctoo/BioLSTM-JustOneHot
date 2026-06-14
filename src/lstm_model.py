"""
The script creates the LSTM model for sequence classification using One-Hot Encoding.
"""
import tensorflow as tf
from tensorflow import keras
from keras.layers import LSTM, Dense, Bidirectional, Dropout, Input
from keras.models import Sequential
from keras.optimizers import Adam
from keras.losses import BinaryCrossentropy

# Parâmetros atualizados para One-Hot e Janela Maior
SEQUENCE_LENGTH = 100
EPOCHS = 10  # Aumentado para compensar a ausência do Embedding
BATCH_SIZE = 100
LSTM_UNITS = 60
LEARNING_RATE = 5e-4
VALIDATION_SPLIT = 0.2

LOSS_FUNCTION = BinaryCrossentropy()
OPTIMIZER = Adam(learning_rate=LEARNING_RATE)

METRICS = [
    'accuracy',
    tf.keras.metrics.Precision(name='precision'),
    tf.keras.metrics.Recall(name='recall'),
    tf.keras.metrics.TruePositives(name='tp'),
    tf.keras.metrics.TrueNegatives(name='tn'),
    tf.keras.metrics.FalsePositives(name='fp'),
    tf.keras.metrics.FalseNegatives(name='fn'),
    tf.keras.metrics.AUC(name='auc')
]

def create_model():
    lstm_model = Sequential([
        # O Input agora recebe a matriz 120x4 do One-Hot
        Input(shape=(SEQUENCE_LENGTH, 4)),

        Bidirectional(LSTM(LSTM_UNITS)),

        # Filtro de decisão padrão
        Dense(32, activation="relu"),
        Dropout(0.3),

        # Veredicto final
        Dense(1, activation="sigmoid")
    ])

    lstm_model.compile(
        optimizer=OPTIMIZER,
        loss=LOSS_FUNCTION,
        metrics=METRICS
    )

    return lstm_model

lstm_model = create_model()
