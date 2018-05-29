import keras
import math
import os
from keras import backend as K

class MyCallback(keras.callbacks.Callback):
        def on_train_begin(self, logs={}):
                current_lr = K.get_value(self.model.optimizer.lr)
                print "Learning rate:", current_lr
                new_lr = current_lr * 0.1
                K.set_value(self.model.optimizer.lr, new_lr)
                new_lr = K.get_value(self.model.optimizer.lr)
                print "New learning rate:", new_lr
                return

class InceptionV4Callback(keras.callbacks.Callback):
        def on_train_begin(self, logs={}):
                current_lr = K.get_value(self.model.optimizer.lr)
                print "Learning rate:", current_lr
                new_lr = current_lr * math.exp(-0.94)
                K.set_value(self.model.optimizer.lr, new_lr)
                new_lr = K.get_value(self.model.optimizer.lr)
                print "New learning rate:", new_lr
                return

class IterationsCallback(keras.callbacks.Callback):

        def __init__(self, validation_generator, validation_steps):
            self.validation_generator = validation_generator
            self.validation_steps = validation_steps
            self.fil = '/scratch/cvn/branch/log/resnet34'

        def on_train_begin(self, logs={}):
                self.losses = []
                self.iteration = 0

                with open(self.fil, 'ar+') as fil:
                    if os.stat(self.fil).st_size == 0:
                        self.losses.append(['iter', 'acc', 'loss', 'val_acc', 'val_loss'])

                    else:
                        self.iteration = int(fil.read().split('\n')[-2].split(' ')[0]) + 1

        def on_batch_end(self, batch, logs={}):
                if self.iteration % 1000 == 0:
                    val_loss, val_acc = self.model.evaluate_generator(self.validation_generator, steps=self.validation_steps)
                    self.losses.append([self.iteration, logs.get('acc'), logs.get('loss'), val_acc, val_loss])
                self.iteration += 1

	def on_epoch_end(self, epoch, logs={}):
                with open(self.fil, 'a') as fil:                    
                    for iteration, acc, loss, val_acc, val_loss in self.losses:
                        fil.write(str(iteration) + ' ' + str(acc) + ' ' + str(loss) + ' ' + str(val_acc) + ' ' + str(val_loss) + '\n')

                self.losses = []
		return


        '''
       	def on_epoch_begin(self, epoch, logs={}):
                current_lr = K.get_value(self.model.optimizer.lr)
                print "Learning rate:", current_lr
                new_lr = 0.02
                K.set_value(self.model.optimizer.lr, new_lr)
                new_lr = K.get_value(self.model.optimizer.lr)
                print "New learning rate:", new_lr
		return

	def on_train_begin(self, logs={}):
                current_lr = K.get_value(self.model.optimizer.lr)
                print "Learning rate:", current_lr
                new_lr = 0.02
                K.set_value(self.model.optimizer.lr, new_lr)
                new_lr = K.get_value(self.model.optimizer.lr)
                print "New learning rate:", new_lr
                return

	def on_train_end(self, logs={}):
		return

	def on_epoch_begin(self, epoch, logs={}):
		return

	def on_epoch_end(self, epoch, logs={}):
		self.losses.append(logs.get('loss'))
		y_pred = self.model.predict(self.validation_data[0])
		self.aucs.append(roc_auc_score(self.validation_data[1], y_pred))
		return

	def on_batch_begin(self, batch, logs={}):
		return

	def on_batch_end(self, batch, logs={}):
		return
        '''
