"""
This is the configuration module.
"""

__version__ = '1.0'
__author__ = 'Saul Alonso-Monsalve'
__email__ = "saul.alonso.monsalve@cern.ch"

import configparser

config = configparser.ConfigParser()

config['random']     = {'seed':'7',
                        'shuffle':'True'}

config['images']     = {'views':'3',
                        'planes':'500',
                        'cells':'500',
                        'standardize':'False',
                        'path':'/scratch/cvn/datasetRaw'}

config['dataset']    = {'uniform':'False',
                        'path':'/scratch/cvn/dataset',
                        'partition_prefix':'/partition',
                        'labels_prefix':'/labels'}

config['log']        = {'path':'/scratch/cvn/log',
                        'prefix':'/train'}

config['model']      = {'architecture':'resnet50',
                        'checkpoint_path':'/scratch/cvn/checkpoint',
                        'checkpoint_prefix':'/model',
                        'checkpoint_save_many':'True',
                        'checkpoint_save_best_only':'False',
                        'checkpoint_period':'1',
                        'parallelize':'False',
                        'gpus':'8',
                        'print_summary':'False',
                        'branches':'False',
                        'outputs':'5'}

config['train']      = {'resume':'False',
                        'weighted_loss_function':'False',
                        'class_weights_prefix':'/class_weights',
                        'lr':'0.1',
                        'momentum':'0.9',
                        'decay':'0.0001',
                        'batch_size':'32',
                        'epochs':'100',
                        'early_stopping_patience':'5',
                        'fraction':'0.98',
                        'max_queue_size':'10'}

config['validation'] = {'batch_size':'32',
                        'fraction':'0.01'}

config['test']       = {'output_path':'/scratch/cvn/output',
                        'output_prefix':'/test_info',
                        'cut_nue':'0.7',
                        'cut_numu':'0.5',
                        'cut_nutau':'0.7',
                        'cut_nc':'0.7',
                        'batch_size':'32',
                        'fraction':'0.01'}

with open('config.ini', 'w') as configfile:
    config.write(configfile)
