import copy
import itertools

import pandas as pd
import pytest

from fgscountrate.fgs_countrate_core import FGS_Countrate
from fgscountrate import utils


def test_compute_fgs_countrate():
    """Test the conversion from GSC magnitudes to FGS countrate"""

    # Test code

    # Test by hand
    id = 'N13I000018'
    fgs = FGS_Countrate(guide_star_id=id, guider=1)
    cr, cr_err, mag, mag_err = fgs.get_fgs_countrate_magnitude()
    assert cr == 1770894.6463908446  # TODO DO THIS BETTER

    # Compare results
    assert 1


def test_fgs_countrate_error():
    """Test the conversion errors due to lacking enough GSC magnitudes"""
    id = 'N13I000018'
    fgs = FGS_Countrate(guide_star_id=id, guider=1)

    # with pytest.raises(NameError) as excinfo:
    #     dataframe = utils.query_gsc(gs_id=id, catalog='GSC241')
    # assert 'No guide stars match' in str(excinfo.value), 'Fake Guide Star ID incorrectly found in catalog'


def test_compute_fgs_magnitude():
    """Test the conversion from GSC magnitudes to FGS magnitude"""

    # Test code

    # Test by hand

    # Compare results
    assert 1
