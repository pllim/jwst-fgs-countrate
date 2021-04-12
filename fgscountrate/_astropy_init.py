# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os

# Create the test function for self test
from astropy.tests.runner import TestRunner
test = TestRunner.make_test_runner_in(os.path.dirname(__file__))
test.__test__ = False
__all__ = ['test']
