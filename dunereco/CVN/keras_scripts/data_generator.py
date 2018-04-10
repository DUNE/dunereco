import numpy as np
import zlib

class DataGenerator(object):
  'Generates data for Keras'

  '''
  Initialization function of the class
  '''

  def __init__(self, cells = 500, planes = 500, views = 3, batch_size = 32, n_labels = 2, interaction_labels = [0,1], 
               neutrino_labels = [], filtered = False, interaction_types = True, images_path = '/', shuffle = True, y_test=[]):
      'Initialization'

      self.cells = cells
      self.planes = planes
      self.views = views
      self.batch_size = batch_size
      self.n_labels = n_labels
      self.interaction_labels = interaction_labels
      self.neutrino_labels = neutrino_labels
      self.filtered = filtered
      self.interaction_types = interaction_types
      self.images_path = images_path
      self.shuffle = shuffle
      self.y_test = y_test
 
  '''
  Goes through the dataset and outputs one batch at a time.
  ''' 

  def generate(self, labels, list_IDs, yield_labels = True):
      'Generates batches of samples'

      # Infinite loop

      while 1:

          # Generate order of exploration of dataset

          indexes = self.__get_exploration_order(list_IDs)

          # Generate batches

          imax = int(len(indexes)/self.batch_size)
          for i in range(imax):

              # Find list of IDs

              list_IDs_temp = [list_IDs[k] for k in indexes[i*self.batch_size:(i+1)*self.batch_size]]

              # Generate data

              if yield_labels:

                  # Train, validation, test

                  X, y = self.__data_generation(labels, list_IDs_temp, yield_labels)

                  yield X, y

              else:

                  # Predictions

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

      X = np.empty((self.batch_size, self.views, self.planes, self.cells))

      if yield_labels:

          # only include the labels when requested (train, validation, and test)

          y = np.empty((self.batch_size), dtype = int)

      # Generate data

      for i, ID in enumerate(list_IDs_temp):

          # Decompress image into pixel NumPy tensor

          if self.filtered:

              # filtered images

              pixels = np.fromstring(zlib.decompress(open(self.images_path + '/' + labels[ID] + '/' + ID + '.txt.gz', 'rb').read()), dtype=np.float64, sep='').reshape(self.views, self.planes, self.cells)
              pixels = pixels.astype('float32')
              pixels/=255

          else:

              # ordinary images

              pixels = np.fromstring(zlib.decompress(open(self.images_path + '/' + labels[ID] + '/' + ID + '.txt.gz', 'rb').read()), dtype=np.uint8, sep='').reshape(self.views, self.planes, self.cells)

          # Store volume

          X[i, :, :, :] = pixels

          # get y value
          
          if self.interaction_types:

              # value from 0 to 13 (12)

              y_value = self.interaction_labels[labels[ID]]

          else:

              # value from 0 to 3

              y_value = self.neutrino_labels[labels[ID]]

          if yield_labels:

              # store class/label (train, validation, and test)

              y[i] = y_value

          else:

              # store actual label (used for the confusion matrix)

              self.y_test.append(y_value)

      if yield_labels:

          # return X and Y (train, validation, and test)

          return X, self.__sparsify(y)

      # return X (predictions)

      return X


  '''
  Please note that Keras only accepts labels written in a binary form 
  (in a 6-label problem, the third label is writtten [0 0 1 0 0 0]), 
  which is why we need the sparsify function to perform this task, 
  should y be a list of numerical values.
  '''

  def __sparsify(self, y):
      'Returns labels in binary NumPy array'

      return np.array([[1 if y[i] == j else 0 for j in range(self.n_labels)] for i in range(y.shape[0])])

