import numpy as np
import zlib, gzip
from string import digits

class DataGenerator(object):
  'Generates data for Keras'

  '''
  Initialization function of the class
  '''

  def __init__(self, cells = 500, planes = 500, views = 3, batch_size = 32, n_labels = 2, interaction_labels = [0,1], neutrino_labels = [], 
               branches = True, standardize = True, filtered = False, interaction_types = True, images_path = '/', shuffle = True, test_values=[]):
      'Initialization'

      self.cells = cells
      self.planes = planes
      self.views = views
      self.batch_size = batch_size
      self.n_labels = n_labels
      self.interaction_labels = interaction_labels
      self.neutrino_labels = neutrino_labels
      self.branches = branches
      self.filtered = filtered
      self.interaction_types = interaction_types
      self.images_path = images_path
      self.standardize = standardize
      self.shuffle = shuffle
      self.test_values = test_values
 
  '''
  Goes through the dataset and outputs one batch at a time.
  ''' 

  def generate(self, labels, list_IDs, yield_labels = True):
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
          
          X = []

          for view in range(self.views):
              X.append(np.empty((self.batch_size, self.planes, self.cells, 1)))

      else:

          # X data should't be a list because there is only one branch

          X = np.empty((self.batch_size, self.planes, self.cells, self.views))
          #X = np.empty((self.batch_size, self.planes, self.cells, 1))

      if yield_labels:

          # only include the labels when requested (train, validation)

          y = np.empty((self.batch_size), dtype = int)

      # Generate data

      for i, ID in enumerate(list_IDs_temp):

          # Decompress image into pixel NumPy tensor

          if self.filtered:

              # filtered images

              with open(self.images_path + '/' + labels[ID] + '/' + ID + '.txt.gz', 'rb') as image_file:      
                  pixels = np.fromstring(zlib.decompress(image_file.read()), dtype=np.float64, sep='').reshape(self.views, self.planes, self.cells)
                  pixels = np.rollaxis(pixels, 0, 3) # from 'channels_first' to 'channels_last'

              #pixels = np.transpose(pixels, (2, 1, 0))

              if self.standardize:

                  # standardize

                  pixels = pixels.astype('float32') # 32-bit precision floating-point pixel image
                  pixels /= 255.                    # pixel range from 0 to 1

          else:

              # ordinary images

              with open(self.images_path + '/' + labels[ID] + '/' + ID + '.txt.gz', 'rb') as image_file:
                  pixels = np.fromstring(zlib.decompress(image_file.read()), dtype=np.uint8, sep='').reshape(self.views, self.planes, self.cells)

              if self.standardize:

                  # standardize the image

                  pixels = pixels.astype('float32') # 32-bit precision floating-point pixel image
                  pixels /= 255.                    # pixel range from 0 to 1

          flavour = str(ID.split('.')[0]).translate(None, digits)          

          # Store volume

          if self.branches:

              for view in range(self.views):

                  X[view][i, :, :, :] = pixels[view, :, :].reshape(self.planes, self.cells, 1)

          else: 

              pixels = np.rollaxis(pixels, 0, 3) # from 'channels_first' to 'channels_last'
              X[i, :, :, :] = pixels

              #X[i, :, :, :] = pixels[2, :, :].reshape(self.planes, self.cells, 1)

          # get y value
          
          if self.interaction_types:

              # value from 0 to 13 (12)

              y_value = self.interaction_labels[labels[ID]]

          else:

              # value from 0 to 3

              y_value = self.neutrino_labels[labels[ID]]

          if flavour[0] == 'a':
              y_value += 13

          if yield_labels:

              # store class/label (train, validation)

              y[i] = y_value

          else:

              # store actual label and energy values (used for the confusion matrix and normalization)

              energy_values = open(self.images_path + '/' + labels[ID] + '/' + ID + '.info', 'rb').readlines()

              self.test_values.append({'flavour': flavour, 
                                       'y_value': y_value,
                                       'fNuEnergy': float(energy_values[0]),
                                       'fRecoNueEnergy': float(energy_values[1]), 
                                       'fRecoNumuEnergy': float(energy_values[2]), 
                                       'fEventWeight': float(energy_values[3])})

      if yield_labels:

          # return X and Y (train, validation)

          return X, self.sparsify3(y)

      # return X (test, predictions)

      return X


  '''
  Please note that Keras only accepts labels written in a binary form 
  (in a 6-label problem, the third label is writtten [0 0 1 0 0 0]), 
  which is why we need the sparsify function to perform this task, 
  should y be a list of numerical values.
  '''

  def sparsify2(self, y):
      'Returns labels in binary NumPy array'

      res1, res2 = [], []

      for i in range(y.shape[0]):

          bi_y1 = [0, 0, 0, 0]  # CC Numu, CC Nue, CC Nutau
          bi_y2 = [0, 0, 0, 0]  # CC QE, CC Res, CC DIS, CC Other

          bi_y1[(y[i] // 4)] = 1 # CC Numu, CC Nue, CC Nutau
          bi_y2[(y[i] %  4)] = 1 # CC QE, CC Res, CC DIS, CC Other

          if y[i] == 12:
              bi_y3 = [-1, -1, -1, -1]

          res1.append(bi_y1)
          res2.append(bi_y2)

      return [np.array(res1), np.array(res2)]

  def sparsify3(self, y):
      'Returns labels in binary NumPy array'

      res1, res2, res3 = [], [], []

      for i in range(y.shape[0]):

          bi_y1 = [0]           # neutrino/antineutrino
          bi_y2 = [0, 0, 0, 0]  # CC Numu, CC Nue, CC Nutau
          bi_y3 = [0, 0, 0, 0]  # CC QE, CC Res, CC DIS, CC Other

          quotient = y[i] // 13

          if quotient > 0:
              y[i] %= 13   # from 0 to 12
              bi_y1[0] = 1 # antineutrino

          bi_y2[(y[i] // 4)] = 1 # CC Numu, CC Nue, CC Nutau
          bi_y3[(y[i] %  4)] = 1 # CC QE, CC Res, CC DIS, CC Other

          if y[i] == 12:
              bi_y1 = [-1]
              bi_y3 = [-1, -1, -1, -1]

          res1.append(bi_y1)
          res2.append(bi_y2)
          res3.append(bi_y3)

      return [np.array(res1), np.array(res2), np.array(res3)]
