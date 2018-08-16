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