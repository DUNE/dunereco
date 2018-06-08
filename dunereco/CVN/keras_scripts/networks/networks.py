import sys
sys.path.append("/home/salonsom/cvn_tensorflow/networks")

from keras.applications.xception import Xception
from keras.applications.vgg16 import VGG16
from keras.applications.vgg19 import VGG19
from keras.applications.resnet50 import ResNet50
from keras.applications.inception_v3 import InceptionV3
from keras.applications.inception_resnet_v2 import InceptionResNetV2
from keras.applications.mobilenet import MobileNet
from keras.applications.densenet import DenseNet121
from keras.applications.densenet import DenseNet169
from keras.applications.densenet import DenseNet201
import my_model, inception_v4, se_resnet, resnet, resnext

from keras.layers import GlobalAveragePooling2D, Dense, Dropout,Activation,Flatten
from keras.layers import Input
from keras.models import Model

'''
################
### Xception ###
################
'''

def _xception(num_classes=13, input_shape=[500, 500, 3], transfer_learning=None):

    print('Architecture: Xception...')

    model = Xception(include_top=True, input_shape=input_shape, weights=None, classes=num_classes)

    return model

'''
################
#### VGG16 #####
################
'''

def _vgg16(num_classes=13, input_shape=[500, 500, 3], transfer_learning=None):

    print('Architecture: VGG16...')

    model = VGG16(include_top=True, input_shape=input_shape, weights=None, classes=num_classes)

    return model

'''
################
#### VGG19 #####
################
'''

def _vgg19(num_classes=13, input_shape=[500, 500, 3], transfer_learning=None):

    print('Architecture: VGG19...')

    model = VGG19(include_top=True, input_shape=input_shape, weights=None, classes=num_classes)

    return model

'''
################
### ResNet18 ###
################
'''

def _resnet18(num_classes=13, input_shape=[500, 500, 3], transfer_learning=None):

    print('Architecture: ResNet18...')

    model = resnet.ResnetBuilder.build_resnet_18(input_shape=input_shape, num_outputs=num_classes)

    return model

'''
################
### ResNet34 ###
################
'''

def _resnet34(num_classes=13, input_shape=[500, 500, 3], transfer_learning=None):

    print('Architecture: ResNet34...')

    model = resnet.ResnetBuilder.build_resnet_34(input_shape=input_shape, num_outputs=num_classes)

    return model

'''
################
### ResNet50 ###
################
'''

def _resnet50(num_classes=13, input_shape=[500, 500, 3], transfer_learning=None):

    print('Architecture: ResNet50...')

    # Ordinary model

    if None == transfer_learning:

        model = ResNet50(include_top=True, input_shape=input_shape, weights=None, classes=num_classes)

    # Transfer learning model

    else:

        pre_model = ResNet50(weights='imagenet', include_top=False)
    
        last_layer = pre_model.output
        x = Flatten()(last_layer)
        x = Dense(num_classes, activation='softmax', name='fc1000')(x)

        model = Model(inputs=pre_model.input, outputs=x)

        if transfer_learning == 'finetuning':
            for layer in model.layers[:-4]:
	        layer.trainable = False
   
    return model

'''
################
## ResNet101 ###
################
'''

def _resnet101(num_classes=13, input_shape=[500, 500, 3], transfer_learning=None):

    print('Architecture: ResNet101...')

    model = resnet.ResnetBuilder.build_resnet_101(input_shape=input_shape, num_outputs=num_classes)

    return model

'''
################
## ResNet152 ###
################
'''

def _resnet152(num_classes=13, input_shape=[500, 500, 3], transfer_learning=None):

    print('Architecture: ResNet152...')

    model = resnet.ResnetBuilder.build_resnet_152(input_shape=input_shape, num_outputs=num_classes)

    return model

'''
################
# InceptionV3 ##
################
'''

def _inceptionv3(num_classes=13, input_shape=[500, 500, 3], transfer_learning=None):

    print('Architecture: InceptionV3...')

    model = InceptionV3(include_top=True, input_shape=input_shape, weights=None, classes=num_classes)

    return model

'''
################
# InceptionV4 ##
################
'''

def _inceptionv4(num_classes=13, input_shape=[500, 500, 3], transfer_learning=None):
    
    print('Architecture: InceptionV4...')

    model = inception_v4.create_model(include_top=True, input_shape=input_shape, weights=None, num_classes=num_classes)

    return model

'''
#####################
# InceptionResNetV2 #
#####################
'''

def _inceptionresnetv2(num_classes=13, input_shape=[500, 500, 3], transfer_learning=None):
    
    print('Architecture: InceptionResNetV2...')

    model = InceptionResNetV2(include_top=True, input_shape=input_shape, weights=None, classes=num_classes)

    return model

'''
################
### ResNeXt ####
################
'''

def _resnext(num_classes=13, input_shape=[500, 500, 3], transfer_learning=None):
    
    print('Architecture: ResNeXt...')

    model = resnext.ResNext(include_top=True, input_shape=input_shape, weights=None, classes=num_classes)

    return model

'''
################
## SEResNet18 ##
################
'''

def _seresnet18(num_classes=13, input_shape=[500, 500, 3], transfer_learning=None):

    print('Architecture: SEResNet18...')

    model = se_resnet.SEResNet18(include_top=True, input_shape=input_shape, weights=None, classes=num_classes)

    return model

'''
################
## SEResNet34 ##
################
'''

def _seresnet34(num_classes=13, input_shape=[500, 500, 3], transfer_learning=None):

    print('Architecture: SEResNet34...')

    model = se_resnet.SEResNet34(include_top=True, input_shape=input_shape, weights=None, classes=num_classes)

    return model

'''
################
## SEResNet50 ##
################
'''

def _seresnet50(num_classes=13, input_shape=[500, 500, 3], transfer_learning=None):

    print('Architecture: SEResNet50...')

    model = se_resnet.SEResNet50(include_top=True, input_shape=input_shape, weights=None, classes=num_classes)

    return model

'''
################
# SEResNet101 ##
################
'''

def _seresnet101(num_classes=13, input_shape=[500, 500, 3], transfer_learning=None):

    print('Architecture: SEResNet101...')

    model = se_resnet.SEResNet101(include_top=True, input_shape=input_shape, weights=None, classes=num_classes)

    return model

'''
################
# SEResNet154 ##
################
'''

def _seresnet154(num_classes=13, input_shape=[500, 500, 3], transfer_learning=None):
    
    print('Architecture: SEResNet154...')

    model = se_resnet.SEResNet154(include_top=True, input_shape=input_shape, weights=None, classes=num_classes)

    return model

'''
################
## MobileNet ###
################
'''

def _mobilenet(num_classes=13, input_shape=[500, 500, 3], transfer_learning=None):
    
    print('Architecture: MobileNet...')

    model = MobileNet(include_top=True, input_shape=input_shape, weights=None, classes=num_classes)

    return model

'''
################
# DenseNet121 ##
################
'''

def _densenet121(num_classes=13, input_shape=[500, 500, 3], transfer_learning=None):

    print('Architecture: DenseNet121...')

    model = DenseNet121(include_top=True, input_shape=input_shape, weights=None, classes=num_classes)

    return model

'''
################
# DenseNet169 ##
################
'''

def _densenet169(num_classes=13, input_shape=[500, 500, 3], transfer_learning=None):

    print('Architecture: DenseNet169...')

    model = DenseNet169(include_top=True, input_shape=input_shape, weights=None, classes=num_classes)

    return model

'''
################
# DenseNet201 ##
################
'''

def _densenet201(num_classes=13, input_shape=[500, 500, 3], transfer_learning=None):

    print('Architecture: DenseNet201...')

    model = DenseNet201(include_top=True, input_shape=input_shape, weights=None, classes=num_classes)

    return model

'''
################
# Custom Model #
################
'''

def _mymodel(num_classes=13, input_shape=[500, 500, 3], transfer_learning=None):

    print('Architecture: Custom model...')

    model = my_model.my_model(input_shape=input_shape, classes=num_classes)

    return model


'''
################
# Create Model #
################
'''

def create_model(network='resnet50', num_classes=13, input_shape=[500, 500, 3], transfer_learning=None):

    network = network.lower()   
 
    if network == 'xception':

        model = _xception(num_classes=num_classes, input_shape=input_shape, transfer_learning=transfer_learning)

    elif network == 'vgg16':

        model = _vgg16(num_classes=num_classes, input_shape=input_shape, transfer_learning=transfer_learning)

    elif network == 'vgg19':

        model = _vgg19(num_classes=num_classes, input_shape=input_shape, transfer_learning=transfer_learning)

    elif network == 'resnet18':

        model = _resnet18(num_classes=num_classes, input_shape=input_shape, transfer_learning=transfer_learning)

    elif network == 'resnet34':

        model = _resnet34(num_classes=num_classes, input_shape=input_shape, transfer_learning=transfer_learning)

    elif network == 'resnet50':

        model = _resnet50(num_classes=num_classes, input_shape=input_shape, transfer_learning=transfer_learning)

    elif network == 'resnet101':

        model = _resnet101(num_classes=num_classes, input_shape=input_shape, transfer_learning=transfer_learning)

    elif network == 'resnet152':

        model = _resnet152(num_classes=num_classes, input_shape=input_shape, transfer_learning=transfer_learning)

    elif network == 'inceptionv3':

        model = _inceptionv3(num_classes=num_classes, input_shape=input_shape, transfer_learning=transfer_learning)

    elif network == 'inceptionv4':

        model = _inceptionv4(num_classes=num_classes, input_shape=input_shape, transfer_learning=transfer_learning)

    elif network == 'inceptionresnetv2':

        model = _inceptionresnetv2(num_classes=num_classes, input_shape=input_shape, transfer_learning=transfer_learning)

    elif network == 'resnext':

        model = _resnext(num_classes=num_classes, input_shape=input_shape, transfer_learning=transfer_learning)

    elif network == 'seresnet18':

        model = _seresnet18(num_classes=num_classes, input_shape=input_shape, transfer_learning=transfer_learning)

    elif network == 'seresnet34':

        model = _seresnet34(num_classes=num_classes, input_shape=input_shape, transfer_learning=transfer_learning)

    elif network == 'seresnet50':

        model = _seresnet50(num_classes=num_classes, input_shape=input_shape, transfer_learning=transfer_learning)

    elif network == 'seresnet101':

        model = _seresnet101(num_classes=num_classes, input_shape=input_shape, transfer_learning=transfer_learning)

    elif network == 'seresnet154':

        model = _seresnet154(num_classes=num_classes, input_shape=input_shape, transfer_learning=transfer_learning)

    elif network == 'mobilenet':

        model = _mobilenet(num_classes=num_classes, input_shape=input_shape, transfer_learning=transfer_learning)

    elif network == 'densenet121':

        model = _densenet121(num_classes=num_classes, input_shape=input_shape, transfer_learning=transfer_learning)

    elif network == 'densenet169':

        model = _densenet169(num_classes=num_classes, input_shape=input_shape, transfer_learning=transfer_learning)

    elif network == 'densenet201':

        model = _densenet201(num_classes=num_classes, input_shape=input_shape, transfer_learning=transfer_learning)

    else:

        model = _mymodel(num_classes=num_classes, input_shape=input_shape, transfer_learning=transfer_learning)

    return model
