import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import sys
import numpy as np
import zlib
import os
from PIL import Image
import glob

# path

for path in glob.iglob('/scratch/datasetRaw/*'):

    index = 0
    print path

    for fil in glob.iglob(path+'/*'):

        if fil[-5:] == ".info":
            continue

        if index >= 20:
            continue

        pixels = np.fromstring(zlib.decompress(open(fil, 'rb').read()), dtype=np.uint8, sep='').reshape(3, 500, 500)
        pixels = np.rollaxis(pixels, 0, 3)

        pixelsaux = np.transpose(pixels, (1,0,2))
        plt.imshow(pixelsaux)
        plt.savefig('/scratch/cvn2/images/numu/event/' + fil.split('/')[-1] + '.png', bbox_inches='tight', pad_inches=0)
        plt.cla()
        index+=1

        for view in range(3):
            p_view = pixels[:, :, view]
            p_view = np.transpose(p_view.reshape(500, 500))

            plt.imshow(p_view)
            plt.savefig('/scratch/cvn2/images/numu/event/' + fil.split('/')[-1] + 'v' + str(view) + '.png', bbox_inches='tight', pad_inches=0)
            plt.cla()
