"""
This is the generator module.
"""

__version__ = '1.0'
__author__ = 'Saul Alonso-Monsalve'
__email__ = "saul.alonso.monsalve@cern.ch"

import numpy as np
import zlib
from string import digits

class DataGenerator(object):

    'Generates data for Keras'

    '''
    Initialization function of the class
    '''
    def __init__(self, cells=500, planes=500, views=3, batch_size=32, branches=True, 
                 outputs=7, standardize=True, images_path = '/', shuffle=True, test_values=[]):
        'Initialization'
        self.cells = cells
        self.planes = planes
        self.views = views
        self.batch_size = batch_size
        self.branches = branches
        self.outputs = outputs
        self.images_path = images_path
        self.standardize = standardize
        self.shuffle = shuffle
        self.test_values = test_values
 
    '''
    Goes through the dataset and outputs one batch at a time.
    ''' 
    def generate(self, labels, list_IDs, yield_labels=True):
        'Generates batches of samples'

        # Infinite loop
        while 1:
            # Generate random order of exploration of dataset (to make each epoch different)
            indexes = self.__get_exploration_order(list_IDs)

            # Generate batches
            imax = int(len(indexes)/self.batch_size) # number of batches

            for i in range(imax):
                 # Find list of IDs for one batch
                 list_IDs_temp = [list_IDs[k] for k in indexes[i*self.batch_size:(i+1)*self.batch_size]]

                 # Generate data
                 if yield_labels:
                     # Train, validation
                     X, y = self.__data_generation(labels, list_IDs_temp, yield_labels)

                     yield X, y
                 else:
                     # Test, predictions
                     X = self.__data_generation(labels, list_IDs_temp, yield_labels)

                     yield X

    '''
    Generates a random order of exploration for a given set of list_IDs. 
    If activated, this feature will shuffle the order in which the examples 
    are fed to the classifier so that batches between epochs do not look alike. 
    Doing so will eventually make our model more robust.
    '''
    def __get_exploration_order(self, list_IDs):
        'Generates order of exploration'

        # Find exploration order
        indexes = np.arange(len(list_IDs))

        if self.shuffle == True:
            np.random.shuffle(indexes)

        return indexes

    '''
    Outputs batches of data and only needs to know about the list of IDs included 
    in batches as well as their corresponding labels.
    '''
    def __data_generation(self, labels, list_IDs_temp, yield_labels):
        'Generates data of batch_size samples' # X : (n_samples, v_size, v_size, v_size, n_channels)

        # Initialization
        if self.branches:
            # X data should be a list of length == branches
            X = [None]*self.views

            for view in range(self.views):
                X[view] = np.empty((self.batch_size, self.planes, self.cells, 1))
        else:
            # X data should't be a list because there is only one branch
            X = np.empty((self.batch_size, self.planes, self.cells, self.views))

        if yield_labels:
            # only include the labels when requested (train, validation)
            if self.outputs == 1:
                y = np.empty((self.batch_size), dtype = int)
            else:
                y = np.empty((self.batch_size, self.outputs), dtype = int)

        # Generate data
        for i, ID in enumerate(list_IDs_temp):
            # Decompress image into pixel NumPy tensor
            with open(self.images_path + '/' + ID.split('.')[0].lstrip('a') + '/images/' + ID + '.gz', 'rb') as image_file:
                pixels = np.fromstring(zlib.decompress(image_file.read()), dtype=np.uint8, sep='').reshape(self.views, self.planes, self.cells)
            #pixels = np.load(self.images_path + '/' + labels[ID] + '/' + ID + '.npy')

            if self.standardize:
                # standardize the image
                pixels = pixels.astype('float32') # 32-bit precision floating-point pixel image
                pixels /= 255.                    # pixel range from 0 to 1

            # Store volume
            if self.branches:
                for view in range(self.views):
                    X[view][i, :, :, :] = pixels[view, :, :].reshape(self.planes, self.cells, 1)
            else: 
                pixels = np.rollaxis(pixels, 0, 3) # from 'channels_first' to 'channels_last'
                X[i, :, :, :] = pixels

            # get y value
            y_value = labels[ID]

            if yield_labels:
                # store class/label (train, validation)
                y[i] = y_value
            else:
                # store actual label and energy values (used for the confusion matrix and normalization)
                with open(self.images_path + '/' + ID.split('.')[0].lstrip('a') + '/info/' + ID + '.info', 'rb') as info_file:
                    energy_values = info_file.readlines()               
                    self.test_values.append({'y_value':y_value,
                                             'fNuEnergy':float(energy_values[1]),
                                             'fLepEnergy':float(energy_values[2]),
                                             'fRecoNueEnergy': float(energy_values[3]), 
                                             'fRecoNumuEnergy': float(energy_values[4]), 
                                             'fEventWeight': float(energy_values[5])})

        if yield_labels:
            # return X and Y (train, validation)
            if self.outputs == 1:
                return X, self.sparsify1(y)
            if self.outputs == 5:
                return X, self.sparsify5(y)
            return X, self.sparsify7(y)

        # return X (test, predictions)
        return X

    '''
    Please note that Keras only accepts labels written in a binary form 
    (in a 6-label problem, the third label is writtten [0 0 1 0 0 0]), 
    which is why we need the sparsify function to perform this task, 
    should y be a list of numerical values.
    '''

    def sparsify1(self, y):
        'Returns labels in binary NumPy array'
        return np.array([[1 if y[i] == j else 1 if y[i]-1 == j and j == 12 else 0 for j in range(13)] for i in range(y.shape[0])])

    def sparsify2(self, y):
        'Returns labels in binary NumPy array'
        res = [None]*2
        res[0] = np.zeros((y.shape[0], 4), dtype=int)
        res[1] = np.zeros((y.shape[0], 4), dtype=int)
 
        for i in range(y.shape[0]):

            res[0][i][(y[i] // 4)] = 1 # CC Numu, CC Nue, CC Nutau
            res[1][i][(y[i] %  4)] = 1 # CC QE, CC Res, CC DIS, CC Other

            if y[i] == 12:
                res[1][i] = [-1, -1, -1, -1]

        return res

    def sparsify3(self, y):
        'Returns labels in binary NumPy array'
        res = [None]*3
        res[0] = np.zeros((y.shape[0], 1), dtype=int)
        res[1] = np.zeros((y.shape[0], 4), dtype=int)
        res[2] = np.zeros((y.shape[0], 4), dtype=int)

        for i in range(y.shape[0]):
            quotient = y[i] // 13

            if quotient > 0:
                y[i] %= 13   # from 0 to 12
                res[0][i][0] = 1 # antineutrino

            res[1][i][(y[i] // 4)] = 1 # CC Numu, CC Nue, CC Nutau
            res[2][i][(y[i] %  4)] = 1 # CC QE, CC Res, CC DIS, CC Other

            if y[i] == 12:
                res[0][i] = [-1]
                res[2][i] = [-1, -1, -1, -1]

        return res

    def normalize(self, value, obj):
        if value == -1 or obj.size == 1:
            obj.fill(value)
        else:
            obj[value] = 1

    def sparsify5(self, y):
        'Returns labels in binary NumPy array'
        res = [None]*self.outputs     

        for i in range(0,len(res)): # flavour, fNProton, fNPion, fNPizero, fNNeutron
            res[i] = np.zeros((y.shape[0], 4), dtype=int)      

        for i in range(y.shape[0]):
            for j in range(len(res)):
                self.normalize(y[i][j], res[j][i])

        return res

    def sparsify7(self, y):
        'Returns labels in binary NumPy array'
        res = [None]*self.outputs
        res[0] = np.zeros((y.shape[0], 1), dtype=int) # fNuPDG

        for i in range(1,len(res)): # flavour, interaction, fNProton, fNPion, fNPizero, fNNeutron
            res[i] = np.zeros((y.shape[0], 4), dtype=int)      

        for i in range(y.shape[0]):
            for j in range(len(res)):
                self.normalize(y[i][j], res[j][i])

        return res
