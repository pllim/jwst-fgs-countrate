# Licensed under a 3-clause BSD style license - see LICENSE.rst

# Packages may add whatever they like to this file, but
# should keep this content at the top.
# ----------------------------------------------------------------------------
from ._astropy_init import *
# ----------------------------------------------------------------------------

import os
import sys
import pkg_resources

module_path = pkg_resources.resource_filename('fgscountrate', '')
setup_path = os.path.normpath(os.path.join(module_path, '../setup.py'))

try:
    with open(setup_path) as f:
        data = f.readlines()

    for line in data:
        if 'VERSION =' in line:
            __version__ = line.split(' ')[-1].replace("'", "").strip()

except FileNotFoundError:
    print('Could not determine fgscountrate version')
    __version__ = '0.0.0'

__minimum_python_version__ = "3.6"


class UnsupportedPythonError(Exception):
    pass

if sys.version_info < tuple((int(val) for val in __minimum_python_version__.split('.'))):
    raise UnsupportedPythonError("fgscountrate does not support Python < {}".format(__minimum_python_version__))

from .fgs_countrate_core import *
from .utils import *
from .conversions import *
