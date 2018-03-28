import configparser

config = configparser.ConfigParser()

config['random']     = {'seed': '7',
                        'shuffle': 'True'}

config['images']     = {'views': '3',
                        'planes': '500',
                        'cells': '500',
                        'labels': {'0':0, '1':1, '2':2, '3':3, '4':4, '5':5, '6':6, '7':7, '8':8, '9':9, '10':10, '11':11, '13':12},
                        'path': '/scratch/devFilesRaw'}

config['dataset']    = {'path': '/scratch/cvn/dataset',
                        'partition_prefix': '/partition',
                        'labels_prefix': '/labels'}

config['log']        = {'path': '/scratch/cvn/log',
                        'prefix': '/log'}

config['model']      = {'checkpoint_path': '/scratch/cvn/checkpoint',
                        'checkpoint_prefix': '/model',
                        'checkpoint_save_many': 'False',
                        'checkpoint_save_best_only': 'True',
                        'checkpoint_period': '1',
                        'print_summary': 'True'}

config['train']      = {'lr': '0.1',
                        'momentum': '0.9',
                        'decay': '0.0000005',
                        'batch_size': '8',
                        'epochs': '10',
                        'early_stopping_patience': '3',
                        'fraction': '0.0006'}

config['validation'] = {'batch_size': '8',
                        'fraction': '0.0002'}

config['test']       = {'batch_size': '8',
                        'fraction': '0.0002'}

with open('config.ini', 'w') as configfile:
    config.write(configfile)
