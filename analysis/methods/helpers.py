import functools
import time
import os
import configparser


def to_str(bytes_or_str):
    if isinstance(bytes_or_str, bytes):
        value = bytes_or_str.decode("utf-8")
    else:
        value = bytes_or_str
    return value


def to_bytes(bytes_or_str):
    if isinstance(bytes_or_str, str):
        value = bytes_or_str.encode('utf-8')
    else:
        value = bytes_or_str
    return value


def all_exist(file_list):
    return all([os.path.isfile(f) for f in file_list])


def process_config(config_file=''):
    """
    by default looks for config file in the same directory as the script
    :param config_file:
    :return:
    """
    if not config_file:
        config_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), "config")
    print(os.path.abspath(__file__))
    config = configparser.ConfigParser()
    config.read(config_file)
    config_dict = {}
    for section in config.sections():
        config_dict[section] = {name: value for name, value in config.items(section)}
    return config_dict


########################################################################################################
# Decorators

def timer(func):
    """ Print the runtime of the decorated funtion
        Quick and dirty, for more precise measuremnt use timeit module
    """
    @functools.wraps(func)
    def wrapper_timer(*args, **kwargs):
        start_time = time.perf_counter()
        value = func(*args, **kwargs)
        end_time = time.perf_counter()
        run_time = end_time - start_time
        print(f"Finished {func.__name__!r} in {run_time:.4f} secs")
        return value
    return wrapper_timer
