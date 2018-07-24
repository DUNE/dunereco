import keras.backend as K
from keras import losses

def masked_loss(y_true, y_pred):
    mask_value = -1
    mask = K.cast(K.not_equal(y_true, mask_value), K.floatx())
        
    return multitask_loss(y_true * mask, y_pred * mask)

def masked_loss_binary(y_true, y_pred):
    mask_value = -1
    mask = K.cast(K.not_equal(y_true, mask_value), K.floatx())
        
    return loss_binary_crossentropy(y_true * mask, y_pred * mask)

def masked_loss_categorical(y_true, y_pred):
    mask_value = -1
    mask = K.cast(K.not_equal(y_true, mask_value), K.floatx())
        
    return loss_categorical_crossentropy(y_true * mask, y_pred * mask)

def loss_binary_crossentropy(y_true, y_pred):
    # Avoid divide by 0
    y_pred = K.clip(y_pred, K.epsilon(), 1 - K.epsilon())

    return losses.binary_crossentropy(y_true, y_pred)

def loss_categorical_crossentropy(y_true, y_pred):
    # Avoid divide by 0
    y_pred = K.clip(y_pred, K.epsilon(), 1 - K.epsilon())

    return losses.categorical_crossentropy(y_true, y_pred)


def multitask_loss(y_true, y_pred):
    # Avoid divide by 0
    y_pred = K.clip(y_pred, K.epsilon(), 1 - K.epsilon())

    # Multi-task loss
    return K.mean(K.sum(- y_true * K.log(y_pred) - (1 - y_true) * K.log(1 - y_pred), axis=1))

def multitask_loss_weighted(y_true, y_pred):
    # Avoid divide by 0
    y_pred = K.clip(y_pred, K.epsilon(), 1 - K.epsilon())

    print "...WEIGHTS..."
    #weights = K.cast([1.5, 1.5, 1.5, 1.5, 0.75, 0.75, 0.75, 0.75], K.floatx())
    #weights = K.cast([1.75, 1.75, 1.75, 1.75, 0.5, 0.5, 0.5, 0.5], K.floatx())
    weights = K.cast([5.0*0.67049114189459802, 5.0*0.6574572346004135, 5.0*2.4884334247030484,  5.0*1.7074017025820696, 
                      0.5*1.1687034551437478,  0.5*0.8630824121263583, 0.5*0.6372003918709082,  0.5*2.4018367377204779], K.floatx())

    # Multi-task loss
    return K.mean(K.sum((- y_true * K.log(y_pred) - (1 - y_true) * K.log(1 - y_pred)) * weights, axis=1))


