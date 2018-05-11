from keras.models import Sequential
from keras.layers import Input, Dense, Dropout, Flatten, BatchNormalization, SeparableConv2D
from keras import regularizers, optimizers
from keras.layers.convolutional import Conv2D, MaxPooling2D, AveragePooling2D

def my_model(input_shape=[500,500,3], classes=3):

    model = Sequential()

    # Convolutional layers

    #model.add(Dropout(0.2, input_shape=input_shape))

    model.add(Conv2D(64, kernel_size=(11,11), strides=4, padding='same', input_shape=input_shape, activation='relu'))
    model.add(MaxPooling2D(pool_size=(4,4), strides=4))
    model.add(BatchNormalization(axis=3))

    model.add(Conv2D(32, kernel_size=(7,7), strides=2, padding='same', activation='relu'))
    model.add(MaxPooling2D(pool_size=(2,2), strides=1))
    model.add(BatchNormalization(axis=3))

    model.add(Dropout(0.2))

    # Flat data to dense layers

    model.add(Flatten())

    # Hiddel layers

    model.add(Dense(1000, 
    #                kernel_regularizer=regularizers.l1_l2(0.001),
    #                activity_regularizer=regularizers.l1_l2(0.001),
                    activation='relu'))

    model.add(Dropout(0.2))

    model.add(Dense(1000,
    #                kernel_regularizer=regularizers.l1_l2(0.001),
    #                activity_regularizer=regularizers.l1_l2(0.001),
                    activation='relu'))

    model.add(Dropout(0.4))

    # Output layer

    model.add(Dense(classes, activation='softmax'))

    return model

