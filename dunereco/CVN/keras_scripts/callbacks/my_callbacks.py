import keras
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
