from pkg_resources import get_distribution

__version__ = get_distribution('g2cpy').version

from .aourw import *
from .genomics import *