import configparser

config = configparser.ConfigParser()

config['random']     = {'seed': '7',
                        'shuffle': 'True'}

config['images']     = {'views': '3',
                        'planes': '500',
                        'cells': '500',
                        'standardize': 'True',
                        'interaction_labels': {'0':0, '1':1, '2':2, '3':3, '4':4, '5':5, '6':6, '7':7, '8':8, '9':9, '10':10, '11':11, '13':12},
                        'neutrino_labels': {'0':1, '1':1, '2':1, '3':1, '4':0, '5':0, '6':0, '7':0, '8':2, '9':2, '10':2, '11':2, '13':3},
                        'path': '/scratch/datasetRaw',
                        'filtered': 'False'}

config['dataset']    = {'uniform': 'False',
                        'interaction_types': 'True', 
                        'path': '/scratch/cvn/dataset',
                        'partition_prefix': '/partition',
                        'labels_prefix': '/labels'}

config['log']        = {'path': '/scratch/cvn/log',
                        'prefix': '/train'}

config['model']      = {'checkpoint_path': '/scratch/cvn/checkpoint',
                        'checkpoint_prefix': '/model',
                        'checkpoint_save_many': 'True',
                        'checkpoint_save_best_only': 'False',
                        'checkpoint_period': '1',
                        'print_summary': 'False'}

config['train']      = {'resume': 'False',
                        'weighted_loss_function': 'False',
                        'class_weights_prefix': '/class_weights',
                        'lr': '0.1',
                        'momentum': '0.9',
                        'decay': '0.0001',
                        'batch_size': '32',
                        'epochs': '100',
                        'early_stopping_patience': '5',
                        'fraction': '0.6'}

config['validation'] = {'batch_size': '32',
                        'fraction': '0.2'}

config['test']       = {'output_path': '/scratch/cvn/output',
                        'output_prefix': '/test_info',
                        'cut_nue': '0.7',
                        'cut_numu': '0.5',
                        'batch_size': '32',
                        'fraction': '0.2'}

with open('config.ini', 'w') as configfile:
    config.write(configfile)
