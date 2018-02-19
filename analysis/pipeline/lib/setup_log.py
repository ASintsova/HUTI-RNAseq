import logging
import logging.config
import configparser
import datetime as dt

def setupLog(analysis_name, logConfig, out_dir):

    log_path = "('{}/{}_{}.log',)".format(out_dir, dt.datetime.now().strftime('%Y_%m_%d_%H_%M_%S'),
                                          analysis_name)
    parser = configparser.ConfigParser()
    parser.read('lib/logConfig')
    parser.set('handler_fileHandler', 'args', log_path)
    parser.set('logger_pipeline', 'qualname', analysis_name)
    with open('lib/logConfig', 'w') as cf:
        parser.write(cf)
    logging.config.fileConfig(logConfig)
    logger = logging.getLogger(analysis_name)
    return logger
