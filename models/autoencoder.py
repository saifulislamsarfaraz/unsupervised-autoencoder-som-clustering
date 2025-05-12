from tensorflow.keras.models import Model
from tensorflow.keras.layers import Input, Dense, Dropout, BatchNormalization
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.regularizers import l2

def build_autoencoder(input_dim, encoding_dim=64):
    input_layer = Input(shape=(input_dim,))
    x = Dense(512, activation='relu', kernel_regularizer=l2(1e-5))(input_layer)
    x = BatchNormalization()(x)
    x = Dropout(0.2)(x)
    encoded = Dense(encoding_dim, activation='relu')(x)
    x = Dense(512, activation='relu')(encoded)
    decoded = Dense(input_dim, activation='linear')(x)

    autoencoder = Model(inputs=input_layer, outputs=decoded)
    encoder = Model(inputs=input_layer, outputs=encoded)
    autoencoder.compile(optimizer=Adam(1e-3), loss='mse')
    return autoencoder, encoder
