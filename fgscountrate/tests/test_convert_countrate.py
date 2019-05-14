import copy
import itertools

import numpy as np
import pandas as pd
import pytest

import fgscountrate
from fgscountrate.fgs_countrate_core import FGSCountrate


def test_compute_countrate_magnitude():
    """Test the conversion from GSC magnitudes to FGS countrate"""

    # Test code

    # Test by hand
    id = 'N13I000018'
    fgs = FGSCountrate(guide_star_id=id, guider=1)
    cr, cr_err, mag, mag_err = fgs.query_fgs_countrate_magnitude()
    assert cr == 1777234.5129574337  # TODO DO THIS BETTER

    # Compare results
    assert 1


def test_errors():
    """Test errors are raised properly """

    gs_id = 'N13I000018'
    guider = 1

    # Test 1: data only includes 2MASS
    fgs = FGSCountrate(guide_star_id=gs_id, guider=guider)
    fgs.gsc_series = fgscountrate.utils.query_gsc(gs_id=id, catalog='GSC241').iloc[0]

    fgs._present_mags = ['tmassJmag', 'tmassHmag', 'tmassKsMag']
    for index in set(fgscountrate.fgs_countrate_core.GSC_BAND_NAMES) - set(fgs._present_mags):
        fgs.gsc_series.loc[index] = -999
    fgs._all_mag_series = fgs.gsc_series.loc[fgscountrate.fgs_countrate_core.GSC_BAND_NAMES]

    with pytest.raises(ValueError) as excinfo:
        fgs.calc_fgs_cr_mag_and_err()
    assert 'Cannot compute' in str(excinfo.value), 'Attempted to compute the FGS countrate & ' \
                                                   'magnitude despite only having the 2MASS bands'

    # Test 2: Guider number is invalid
    guider = 3
    fgs = FGSCountrate(guide_star_id=gs_id, guider=guider)
    fgs.gsc_series = fgscountrate.utils.query_gsc(gs_id=id, catalog='GSC241').iloc[0]

    with pytest.raises(ValueError) as excinfo:
        fgs.calc_fgs_cr_mag_and_err()
    assert '1 or 2' in str(excinfo.value), 'Allowed invalid guider number to pass'


def test_output_options():
    """
    Test the output options for calc_fgs_cr_mag_and_err()
        and _calc_fgs_cr_mag() are as expected
    """

    gs_id = 'N13I000018'
    guider = 2
    fgs = FGSCountrate(guide_star_id=gs_id, guider=guider)
    fgs.gsc_series = fgscountrate.utils.query_gsc(gs_id=gs_id, catalog='GSC241').iloc[0]
    fgs._present_mags = fgscountrate.fgs_countrate_core.GSC_BAND_NAMES
    fgs._all_mag_series = fgs.gsc_series.loc[fgs._present_mags]
    mag_err_list = [fgs.gsc_series[ind + 'Err'] for ind in fgs._all_mag_series.index]
    fgs._all_mag_err_series = pd.Series(mag_err_list, index=fgs._all_mag_series.index)

    # Test output from calc_fgs_cr_mag_and_err()
    return_list = fgs.calc_fgs_cr_mag_and_err()
    assert len(return_list) == 4

    # Test output from _calc_fgs_cr_mag()
    band_series = fgs._all_mag_series
    guider_throughput = fgscountrate.fgs_countrate_core.THROUGHPUT_G2
    guider_gain = fgscountrate.fgs_countrate_core.CR_CONVERSION_G2

    # Case 1: Only Countrate
    return_list = fgs._calc_fgs_cr_mag(to_compute='countrate',
                                       band_series=band_series, guider_throughput=guider_throughput,
                                       guider_gain=guider_gain, return_dataframe=False)
    assert len(return_list) == 1
    assert pytest.approx(return_list[0], fgs.fgs_countrate, 1e-5)

    # Case 2: Only Magnitude
    return_list = fgs._calc_fgs_cr_mag(to_compute='magnitude',
                                       band_series=band_series, guider_throughput=guider_throughput,
                                       guider_gain=guider_gain, return_dataframe=False)
    assert len(return_list) == 1
    assert pytest.approx(return_list[0], fgs.fgs_magnitude, 1e-5)

    # Case 3: Both
    return_list = fgs._calc_fgs_cr_mag(to_compute='both',
                                       band_series=band_series, guider_throughput=guider_throughput,
                                       guider_gain=guider_gain, return_dataframe=False)
    assert len(return_list) == 2
    assert pytest.approx(return_list[0], fgs.fgs_countrate, 1e-5)
    assert pytest.approx(return_list[1], fgs.fgs_magnitude, 1e-5)

    # Case 4: Countrate + Dataframe Only
    return_list = fgs._calc_fgs_cr_mag(to_compute='countrate',
                                       band_series=band_series, guider_throughput=guider_throughput,
                                       guider_gain=guider_gain, return_dataframe=True)
    assert len(return_list) == 2
    assert pytest.approx(return_list[0], fgs.fgs_countrate, 1e-5)
    np.testing.assert_array_almost_equal(return_list[1].values.flatten().tolist(),
                                         fgs.band_dataframe.values.flatten().tolist(), 5)

    # Case 5: Magnitude + Dataframe Only
    return_list = fgs._calc_fgs_cr_mag(to_compute='magnitude',
                                       band_series=band_series, guider_throughput=guider_throughput,
                                       guider_gain=guider_gain, return_dataframe=True)
    assert len(return_list) == 2
    assert pytest.approx(return_list[0], fgs.fgs_magnitude, 1e-5)
    np.testing.assert_array_almost_equal(return_list[1].values.flatten().tolist(),
                                         fgs.band_dataframe.values.flatten().tolist(), 5)

    # Case 6: Countrate, Magnitude, and Dataframe
    return_list = fgs._calc_fgs_cr_mag(to_compute='both',
                                       band_series=band_series, guider_throughput=guider_throughput,
                                       guider_gain=guider_gain, return_dataframe=True)
    assert len(return_list) == 3
    assert pytest.approx(return_list[0], fgs.fgs_countrate, 1e-5)
    assert pytest.approx(return_list[1], fgs.fgs_magnitude, 1e-5)
    np.testing.assert_array_almost_equal(return_list[2].values.flatten().tolist(),
                                         fgs.band_dataframe.values.flatten().tolist(), 5)
