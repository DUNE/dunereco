import matplotlib
matplotlib.use('Agg')

from keras.models import Model, load_model
import matplotlib.pyplot as plt
import sys
import numpy as np
import zlib
import os
from PIL import Image

# path

IMAGES_PATH = '/scratch/cvn2/images/nc'
MODEL_PATH  = '/scratch/cvn2/checkpoint/modelaux-09-0.77.h5'
IMAGE_PATH  = '/scratch/datasetRaw/9/nue2.list8737.info'

if not os.path.exists(IMAGES_PATH + '/event'):
    os.makedirs(IMAGES_PATH + '/event')
if not os.path.exists(IMAGES_PATH + '/feature_maps'):
    os.makedirs(IMAGES_PATH + '/feature_maps')
if not os.path.exists(IMAGES_PATH + '/filters'):
    os.makedirs(IMAGES_PATH + '/filters')

# load model

model = load_model(MODEL_PATH)

# print model

#print model.summary()
#print model2.summary()

# plot the filters

conv_num = 0

pixels = np.fromstring(zlib.decompress(open(IMAGE_PATH, 'rb').read()), dtype=np.uint8, sep='').reshape(3, 500, 500)
pixels = np.rollaxis(pixels, 0, 3)

# save event

pixelsaux = np.transpose(pixels, (1,0,2))
plt.imshow(pixelsaux)
plt.axis('off')
plt.savefig(IMAGES_PATH + '/event/event.pdf', bbox_inches='tight', pad_inches=0)
plt.cla()

for view in range(3):
    p_view = pixels[:, :, view]
    p_view = np.transpose(p_view.reshape(500, 500))

    plt.imshow(p_view)
    plt.axis('off')
    plt.savefig(IMAGES_PATH + '/event/view' + str(view) + '.pdf', bbox_inches='tight', pad_inches=0)
    plt.cla()

X = pixels.reshape(1, 500, 500, 3)

index = 0

for layer in model.layers:

    print "Layer name:", layer.name

    if layer.name[:5] == "input":
        print "    Input shape:", model.input.shape

    if layer.name[:4] != "conv":
        continue
   
    print "    Shape of the filters:", layer.get_weights()[0].shape

    weights = layer.get_weights()[0]

    '''
    # filters

    print "        saving filters..."

    for filt in range(weights.shape[3]):
        for channel in range(weights.shape[2]):
            kernel = weights[:,:,channel,filt]

            kernel = np.transpose(kernel.reshape(kernel.shape[0], kernel.shape[1]))

            if(kernel.min() < 0):
                kernel = kernel - kernel.min()

            kernel /= kernel.max()

            if not os.path.exists(IMAGES_PATH + '/filters/' + str(index)):
                os.makedirs(IMAGES_PATH + '/filters/' + str(index))

            plt.imshow(kernel)
            plt.savefig(IMAGES_PATH + '/filters/' + str(index) + '/' + str(filt) + '_' + str(channel)  + '.pdf', bbox_inches='tight', pad_inches=0)
            plt.cla()
    '''

    model2 = Model(inputs=model.input, outputs=model.get_layer(layer.name).output)
    feature_maps = model2.predict(X)

    print "    Shape of the feature maps:", feature_maps.shape

    # feature maps

    print "        saving feature maps..."

    for i in range(feature_maps.shape[3]):

        feature_map = feature_maps[:,:,:,i]
        feature_map = np.transpose(feature_map.reshape(feature_maps.shape[1], feature_maps.shape[2]))
        feature_map = np.asarray(feature_map, dtype=np.uint8)

        if not os.path.exists(IMAGES_PATH + '/feature_maps/' + str(index)):
            os.makedirs(IMAGES_PATH + '/feature_maps/' + str(index))

        plt.imshow(feature_map)
        plt.axis('off')
        plt.savefig(IMAGES_PATH + '/feature_maps/' + str(index) + '/fm' + str(i) + '.pdf', bbox_inches='tight', pad_inches=0)
        plt.cla()

    index += 1
