import keras
import math
import os
import numpy as np
from keras import backend as K

class MultiGPUCheckpointCallback(keras.callbacks.Callback):

    def __init__(self, filepath, base_model, monitor='val_loss', verbose=0,
                 save_best_only=False, save_weights_only=False,
                 mode='auto', period=1):
        super(MultiGPUCheckpointCallback, self).__init__()
        self.base_model = base_model
        self.monitor = monitor
        self.verbose = verbose
        self.filepath = filepath
        self.save_best_only = save_best_only
        self.save_weights_only = save_weights_only
        self.period = period
        self.epochs_since_last_save = 0

        if mode not in ['auto', 'min', 'max']:
            warnings.warn('ModelCheckpoint mode %s is unknown, '
                          'fallback to auto mode.' % (mode),
                          RuntimeWarning)
            mode = 'auto'

        if mode == 'min':
            self.monitor_op = np.less
            self.best = np.Inf
        elif mode == 'max':
            self.monitor_op = np.greater
            self.best = -np.Inf
        else:
            if 'acc' in self.monitor or self.monitor.startswith('fmeasure'):
                self.monitor_op = np.greater
                self.best = -np.Inf
            else:
                self.monitor_op = np.less
                self.best = np.Inf

    def on_epoch_end(self, epoch, logs=None):

        logs = logs or {}
        self.epochs_since_last_save += 1
        if self.epochs_since_last_save >= self.period:
            self.epochs_since_last_save = 0
            filepath = self.filepath.format(epoch=epoch + 1, **logs)
            if self.save_best_only:
                current = logs.get(self.monitor)
                if current is None:
                    warnings.warn('Can save best model only with %s available, '
                                  'skipping.' % (self.monitor), RuntimeWarning)
                else:
                    if self.monitor_op(current, self.best):
                        if self.verbose > 0:
                            print('Epoch %05d: %s improved from %0.5f to %0.5f,'
                                  ' saving model to %s'
                                  % (epoch + 1, self.monitor, self.best,
                                     current, filepath))
                        self.best = current
                        if self.save_weights_only:
                            self.base_model.save_weights(filepath, overwrite=True)
                        else:
                            self.base_model.save(filepath, overwrite=True)
                    else:
                        if self.verbose > 0:
                            print('Epoch %05d: %s did not improve' %
                                  (epoch + 1, self.monitor))
            else:
                if self.verbose > 0:
                    print('Epoch %05d: saving model to %s' % (epoch + 1, filepath))
                if self.save_weights_only:
                    self.base_model.save_weights(filepath, overwrite=True)
                else:
                    self.base_model.save(filepath, overwrite=True)

def detachmodel(m):
    """ Detach model trained on GPUs from its encapsulation
    # Arguments
        :param m: obj, keras model
    # Returns
        :return: obj, keras model
    """

    for l in m.layers:
        if l.name == 'resnext':
            return l

    return m

class ModelCheckpointDetached(keras.callbacks.Callback):

    """ Save detached from multi-GPU encapsulation model
    (very small) modification from https://github.com/fchollet/keras/blob/master/keras/callbacks.py#L331

    `filepath` can contain named formatting options,
    which will be filled the value of `epoch` and
    keys in `logs` (passed in `on_epoch_end`).

    For example: if `filepath` is `weights.{epoch:02d}-{val_loss:.2f}.hdf5`,
    then the model checkpoints will be saved with the epoch number and
    the validation loss in the filename.

    # Arguments
        filepath: string, path to save the model file.
        monitor: quantity to monitor.
        verbose: verbosity mode, 0 or 1.
        save_best_only: if `save_best_only=True`,
            the latest best model according to
            the quantity monitored will not be overwritten.
        mode: one of {auto, min, max}.
            If `save_best_only=True`, the decision
            to overwrite the current save file is made
            based on either the maximization or the
            minimization of the monitored quantity. For `val_acc`,
            this should be `max`, for `val_loss` this should
            be `min`, etc. In `auto` mode, the direction is
            automatically inferred from the name of the monitored quantity.
        save_weights_only: if True, then only the model's weights will be
            saved (`model.save_weights(filepath)`), else the full model
            is saved (`model.save(filepath)`).
        period: Interval (number of epochs) between checkpoints.
    """

    def __init__(self, filepath, monitor='val_loss', verbose=0,
                 save_best_only=False, save_weights_only=False,
                 mode='auto', period=1):
        super(ModelCheckpointDetached, self).__init__()
        self.monitor = monitor
        self.verbose = verbose
        self.filepath = filepath
        self.save_best_only = save_best_only
        self.save_weights_only = save_weights_only
        self.period = period
        self.epochs_since_last_save = 0

        if mode not in ['auto', 'min', 'max']:
            warnings.warn('ModelCheckpoint mode %s is unknown, '
                          'fallback to auto mode.' % mode, RuntimeWarning)
            mode = 'auto'

        if mode == 'min':
            self.monitor_op = np.less
            self.best = np.Inf
        elif mode == 'max':
            self.monitor_op = np.greater
            self.best = -np.Inf
        else:
            if 'acc' in self.monitor or self.monitor.startswith('fmeasure'):
                self.monitor_op = np.greater
                self.best = -np.Inf
            else:
                self.monitor_op = np.less
                self.best = np.Inf

    def on_epoch_end(self, epoch, logs=None):
        logs = logs or {}
        self.epochs_since_last_save += 1
        if self.epochs_since_last_save >= self.period:
            self.epochs_since_last_save = 0
            filepath = self.filepath.format(epoch=epoch, **logs)
            if self.save_best_only:
                current = logs.get(self.monitor)
                if current is None:
                    warnings.warn('Can save best model only with %s available, '
                                  'skipping.' % self.monitor, RuntimeWarning)
                else:
                    if self.monitor_op(current, self.best):
                        if self.verbose > 0:
                            print('Epoch %05d: %s improved from %0.5f to %0.5f,'
                                  ' saving model to %s'
                                  % (epoch, self.monitor, self.best,
                                     current, filepath))
                        self.best = current
                        if self.save_weights_only:
                            detachmodel(self.model).save_weights(filepath, overwrite=True)
                        else:
                            detachmodel(self.model).save(filepath, overwrite=True)
                    else:
                        if self.verbose > 0:
                            print('Epoch %05d: %s did not improve' %
                                  (epoch, self.monitor))
            else:
                if self.verbose > 0:
                    print('Epoch %05d: saving model to %s' % (epoch, filepath))
                if self.save_weights_only:
                    detachmodel(self.model).save_weights(filepath, overwrite=True)
                else:
                    detachmodel(self.model).save(filepath, overwrite=True)

class MyCallback(keras.callbacks.Callback):
        def on_epoch_end(self, epoch, logs={}):
                current_lr = K.get_value(self.model.optimizer.lr)
                print "Learning rate:", current_lr
                new_lr = current_lr * 0.95
                K.set_value(self.model.optimizer.lr, new_lr)
                new_lr = K.get_value(self.model.optimizer.lr)
                print "New learning rate:", new_lr
                return

class InceptionV4Callback(keras.callbacks.Callback):
        def on_train_begin(self, logs={}):
                current_lr = K.get_value(self.model.optimizer.lr)
                print "Learning rate:", current_lr
                new_lr = current_lr * 0.94
                K.set_value(self.model.optimizer.lr, new_lr)
                new_lr = K.get_value(self.model.optimizer.lr)
                print "New learning rate:", new_lr
                return

class InceptionV3Callback(keras.callbacks.Callback):
	def on_epoch_end(self, epoch, logs={}):
                if epoch % 2 == 1:
                    current_lr = K.get_value(self.model.optimizer.lr)
                    print "Learning rate:", current_lr
                    new_lr = current_lr * 0.94
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
