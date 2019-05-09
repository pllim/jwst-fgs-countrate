import copy
import itertools

import pandas as pd
import pytest

from fgscountrate.fgs_countrate_core import FGS_Countrate


def test_compute_fgs_countrate():
    """Test the conversion from GSC magnitudes to FGS countrate"""

    # Test code

    # Test by hand
    id = 'N13I000018'
    fgs = FGS_Countrate(guide_star_id=id, guider=1)
    cr = fgs.get_fgs_countrate()
    assert cr == 846062065.392

    # Compare results
    assert 1


def test_fgs_countrate_error():
    """Test the conversion errors due to lacking enough GSC magnitudes"""
    id = 'N13I000018'
    fgs = FGS_Countrate(guide_star_id=id, guider=1)

    with pytest.raises(NameError) as excinfo:
        dataframe = utils.query_gsc(gs_id=id, catalog='GSC241')
    assert 'No guide stars match' in str(excinfo.value), 'Fake Guide Star ID incorrectly found in catalog'


def test_compute_fgs_magnitude():
    """Test the conversion from GSC magnitudes to FGS magnitude"""

    # Test code

    # Test by hand

    # Compare results
    assert 1
