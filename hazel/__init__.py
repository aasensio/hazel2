from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("hazel2")
except PackageNotFoundError:
    __version__ = "unknown"
__author__ = "Andres Asensio Ramos"

from .atmosphere import *
from .model import *
from .configuration import *
from .spectrum import *
from .multiprocess import *
from .tools import *
from .io_hazel import *
from .sir import *
from .exceptions import *
from .nlte_weights import *
from .util import *
from .analysis import *
from . import codes

# try:
#     # import torch
#     # import torch_geometric
#     from .graphnet import *
#     from .forward_nn import *    
# except:
#     pass

try:
    from .forward_nn_transformer import *
except:
    pass
