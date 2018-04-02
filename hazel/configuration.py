import numpy as np
import matplotlib.pyplot as pl
from configobj import ConfigObj
from ipdb import set_trace as stop
__all__ = ['Configuration']

def _lower_to_sep(string, separator='='):
    line=string.partition(separator)
    string=str(line[0]).lower()+str(line[1])+str(line[2])
    return string

class Configuration(object):

    def __init__(self, filename):

        f = open(filename, 'r')
        tmp = f.readlines()
        f.close()

        input_lower = ['']

        for l in tmp:
            input_lower.append(_lower_to_sep(l)) # Convert keys to lowercase

        # Parse configuration file
        self.config_dict = ConfigObj(input_lower)