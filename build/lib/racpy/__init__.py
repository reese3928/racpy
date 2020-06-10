"""Top-level package for racpy."""

__author__ = """Xu Ren"""
__email__ = 'xuren2120@gmail.com'
__version__ = '0.1.4'

from .core import RNAAgeCalc
from .makeplot import makeplot

import pandas as pd
import os

location = os.path.dirname(os.path.realpath(__file__))
fpkm = pd.read_csv(os.path.join(location, "internal_data", "fpkm.csv"),
                   index_col=0)
rawcount = pd.read_csv(os.path.join(location, "internal_data", "rawcount.csv"),
                       index_col=0)
