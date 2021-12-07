import copy
import itertools
import warnings

import numpy as np
import pandas as pd
import pytest

import fgscountrate
from fgscountrate.fgs_countrate_core import FGSCountrate
from fgscountrate.constants import GSC_BAND_NAMES

values = ['N13I000018 ', 420900912, 273.206729760604, 65.5335149359777,
          8.3030233068735e-05, 0.000185964552890292, 14.9447, 0.285722,
          14.0877, 0.29279299999999997, 13.7468, 0.239294, 13.33899974823,
          0.025000000372529, 12.9930000305176, 0.0270000007003546,
          12.9010000228882, 0.0270000007003546, 15.78594, 0.005142466,
          14.654670000000001, 0.003211281, 14.27808, 0.0032733809999999997,
          14.14432, 0.003414216, 14.106670000000001, 0.00433389]
index = ['hstID', 'gsc1ID', 'ra', 'dec', 'raErr', 'decErr',
         'JpgMag', 'JpgMagErr', 'FpgMag', 'FpgMagErr', 'NpgMag', 'NpgMagErr',
         'tmassJMag', 'tmassJMagErr', 'tmassHMag', 'tmassHMagErr', 'tmassKsMag', 'tmassKsMagErr',
         'SDSSuMag', 'SDSSuMagErr', 'SDSSgMag', 'SDSSgMagErr', 'SDSSrMag', 'SDSSrMagErr',
         'SDSSiMag', 'SDSSiMagErr', 'SDSSzMag', 'SDSSzMagErr']
GSC_SERIES = pd.Series(values, index=index)


test_query_fgs_countrate_magnitude_parameters = ['GSC242', 'GSC241']
@pytest.mark.parametrize('gsc', test_query_fgs_countrate_magnitude_parameters)
def test_query_fgs_countrate_magnitude(gsc):
    """
    Test this function runs smoothly and doesn't error
    for multiple guide star catalogs
    (not tested anywhere else)
    """

    # A case with all bands/uncertainties present
    gs_id = 'N13I000018'
    guider = 1
    fgs = FGSCountrate(guide_star_id=gs_id, guider=guider)
    cr, cr_err, mag, mag_err = fgs.query_fgs_countrate_magnitude(catalog=gsc)
    assert cr > 0
    assert mag > 0
    if any(i == -999 for i in fgs._all_queried_mag_series.values):
        warnings.warn('GS ID N13I000018 no longer behaves as originally expected. Test must be updated')

    # A case with missing all 3 tmass bands
    gs_id = 'N94D006388'
    guider = 1
    fgs = FGSCountrate(guide_star_id=gs_id, guider=guider)
    cr, cr_err, mag, mag_err = fgs.query_fgs_countrate_magnitude()
    assert any(fgs.band_dataframe.at[i, 'Signal'] > 0 for i in ['tmassJMag', 'tmassHMag', 'tmassKsMag'])
    assert cr > 0
    assert mag > 0
    if any(i != -999 for i in fgs.gsc_series[['tmassJMag', 'tmassHMag', 'tmassKsMag']].values):
        warnings.warn('GS ID N94D006388 no longer behaves as originally expected. Test must be updated')

    # A case with missing uncertainty data
    gs_id = 'N13I018276'
    guider = 1
    fgs = FGSCountrate(guide_star_id=gs_id, guider=guider)
    cr, cr_err, mag, mag_err = fgs.query_fgs_countrate_magnitude()
    assert fgs.k_mag_err > 0  # this value should get reset
    assert cr > 0
    assert mag > 0
    if fgs.gsc_series['tmassKsMagErr'] != -999:
        warnings.warn('GS ID N13I018276 no longer behaves as originally expected. Test must be updated')


def test_compute_countrate_magnitude():
    """
    Test the conversion from GSC magnitudes to FGS countrate
    returns values as expected
    """

    gs_id = 'N13I000018'
    guider = 1
    fgs = FGSCountrate(guide_star_id=gs_id, guider=guider)

    # Reset data to a set of constant, fake data
    fgs.gsc_series = copy.copy(GSC_SERIES)

    # Convert to JHK magnitudes
    fgs.j_mag, fgs.j_mag_err, fgs.h_mag, fgs.h_mag_err, fgs.k_mag, fgs.k_mag_err = \
        fgs.calc_jhk_mag(fgs.gsc_series)

    # Compute FGS countrate and magnitude
    cr, cr_err, mag, mag_err = fgs.calc_fgs_cr_mag_and_err()

    assert np.isclose(cr, 1786779.2896366853, 1e-5)
    assert np.isclose(cr_err, 153161.72059228245, 1e-5)
    assert np.isclose(mag, 13.304318358662279, 1e-5)
    assert np.isclose(mag_err, 0.665287435671057, 1e-5)


def test_convert_cr_to_fgs_mag():
    """Test count rate to magnitude conversion helper function"""
    # Numbers come from case with all bands
    countrate = 1777234.5129574337
    expected_mag = 13.310964314752303
    mag = fgscountrate.convert_cr_to_fgs_mag(countrate, guider=1)
    assert np.isclose(mag, expected_mag, 1e-5)

    # Numbers come from case with missing bands
    countrate = 1815659.5085523769
    expected_mag = 13.28774013985303
    mag = fgscountrate.convert_cr_to_fgs_mag(countrate, guider=1)
    assert np.isclose(mag, expected_mag, 1e-5)


def test_convert_fgs_mag_to_cr():
    """Test magnitude to count rate conversion helper function"""
    # Numbers come from case with all bands
    magnitude = 13.310964314752303
    expected_cr = 1777234.5129574337
    cr = fgscountrate.convert_fgs_mag_to_cr(magnitude, guider=1)
    assert np.isclose(cr, expected_cr, 1)

    # Numbers come from case with missing bands
    magnitude = 13.28774013985303
    expected_cr = 1815659.5085523769
    cr = fgscountrate.convert_fgs_mag_to_cr(magnitude, guider=1)
    assert np.isclose(cr, expected_cr, 1)


def test_band_missing():
    """Test that when a band (SDSS_g) is missing, it's signal is set to 0"""

    gs_id = 'N13I000018'
    guider = 1
    fgs = FGSCountrate(guide_star_id=gs_id, guider=guider)

    # Reset data to a set of constant, fake data with SDSS_g missing
    fgs.gsc_series = copy.copy(GSC_SERIES)
    fgs.gsc_series['SDSSgMag'] = -999
    fgs.gsc_series['SDSSgMagErr'] = -999

    # Convert to JHK magnitudes
    fgs.j_mag, fgs.j_mag_err, fgs.h_mag, fgs.h_mag_err, fgs.k_mag, fgs.k_mag_err = \
        fgs.calc_jhk_mag(fgs.gsc_series)

    # Compute FGS countrate and magnitude to get fgs.band_dataframe attribute
    _ = fgs.calc_fgs_cr_mag_and_err()

    # Check Mag, ABMag, Flux, and Signal = -999
    assert fgs.survey == 'sdss'
    assert fgs.band_dataframe.at['SDSSgMag', 'Mag'] == -999
    assert fgs.band_dataframe.at['SDSSgMag', 'ABMag'] == -999
    assert fgs.band_dataframe.at['SDSSgMag', 'Flux'] == -999
    assert fgs.band_dataframe.at['SDSSgMag', 'Signal'] == -999


def test_dim_limits_partialsdss():
    """Test that when some SDSS bands are below the dim limits, those bands
    are not included in the JHK or count rate/magnitude calculations, but
    the remaining SDSS bands are used."""

    gs_id = 'N13I000018'
    guider = 1
    fgs = FGSCountrate(guide_star_id=gs_id, guider=guider)

    # Reset data to a set of constant, fake data with TMASS not present and 2 SDSS bands set to dim
    fgs.gsc_series = copy.copy(GSC_SERIES)
    fgs.gsc_series['tmassJMag'] = -999
    fgs.gsc_series['tmassHMag'] = -999
    fgs.gsc_series['tmassKsMag'] = -999
    fgs.gsc_series['SDSSgMag'] = 24  # dim
    fgs.gsc_series['SDSSrMag'] = 24  # dim
    fgs.gsc_series['SDSSzMag'] = 17
    fgs.gsc_series['SDSSiMag'] = 17

    # Convert to JHK magnitudes
    fgs.j_mag, fgs.j_mag_err, fgs.h_mag, fgs.h_mag_err, fgs.k_mag, fgs.k_mag_err = \
        fgs.calc_jhk_mag(fgs.gsc_series)

    # Check conversion method
    assert fgs.j_convert_method == 'convert_sdssiz_to_jhk'
    assert fgs.h_convert_method == 'convert_sdssiz_to_jhk'
    assert fgs.k_convert_method == 'convert_sdssiz_to_jhk'

    # Compute FGS countrate and magnitude to get fgs.band_dataframe attribute
    _ = fgs.calc_fgs_cr_mag_and_err()

    # Check Mag, ABMag, Flux, and Signal = -999 for g and r
    assert fgs.survey == 'sdss'

    assert fgs.band_dataframe.at['SDSSgMag', 'ABMag'] == -999
    assert fgs.band_dataframe.at['SDSSgMag', 'Flux'] == -999
    assert fgs.band_dataframe.at['SDSSgMag', 'Signal'] == -999

    assert fgs.band_dataframe.at['SDSSrMag', 'ABMag'] == -999
    assert fgs.band_dataframe.at['SDSSrMag', 'Flux'] == -999
    assert fgs.band_dataframe.at['SDSSrMag', 'Signal'] == -999


def test_dim_limits_allsdss():
    """Test that when all SDSS bands are below the dim limits, they are not included in the
    JHK or count rate/magnitude calculations, and GSC2 bands are used instead."""

    gs_id = 'N13I000018'
    guider = 1
    fgs = FGSCountrate(guide_star_id=gs_id, guider=guider)

    # Reset data to a set of constant, fake data with TMASS not present and all SDSS bands set to dim
    fgs.gsc_series = copy.copy(GSC_SERIES)
    fgs.gsc_series['tmassJMag'] = -999
    fgs.gsc_series['tmassHMag'] = -999
    fgs.gsc_series['tmassKsMag'] = -999
    fgs.gsc_series['SDSSgMag'] = 24
    fgs.gsc_series['SDSSrMag'] = 24
    fgs.gsc_series['SDSSzMag'] = 24
    fgs.gsc_series['SDSSiMag'] = 24

    # Convert to JHK magnitudes
    fgs.j_mag, fgs.j_mag_err, fgs.h_mag, fgs.h_mag_err, fgs.k_mag, fgs.k_mag_err = \
        fgs.calc_jhk_mag(fgs.gsc_series)

    # Check conversion method
    assert fgs.j_convert_method == 'convert_gsc2bjin_to_jhk'
    assert fgs.h_convert_method == 'convert_gsc2bjin_to_jhk'
    assert fgs.k_convert_method == 'convert_gsc2bjin_to_jhk'

    # Check that no SDSS bands made it into the present bands list
    assert ['sdss' not in substring.lower() for substring in fgs._present_queried_mags]
    assert ['sdss' not in substring.lower() for substring in fgs._present_calculated_mags]

    # Compute FGS countrate and magnitude to get fgs.band_dataframe attribute
    _ = fgs.calc_fgs_cr_mag_and_err()

    # Check this calculation should only be done with GSC bands
    assert fgs.survey == 'gsc2'
    assert fgs.band_dataframe.at['JpgMag', 'ABMag'] != -999
    assert fgs.band_dataframe.at['FpgMag', 'ABMag'] != -999
    assert fgs.band_dataframe.at['NpgMag', 'ABMag'] != -999


def test_sdss_or_gsc():
    """Test that only SDSS or GSC data is used to compute CR and Mag, and not a combination of the two"""

    gs_id = 'N13I000018'
    guider = 1

    # Test SDSS only - set tmass to missing (sdss will be chosen over gsc already)
    fgs = FGSCountrate(guide_star_id=gs_id, guider=guider)
    fgs.gsc_series = copy.copy(GSC_SERIES)
    fgs.gsc_series['tmassJMag'] = -999
    fgs.gsc_series['tmassHMag'] = -999
    fgs.gsc_series['tmassKsMag'] = -999
    fgs.calc_jhk_mag(fgs.gsc_series)
    _ = fgs.calc_fgs_cr_mag_and_err()
    # check gsc bands have been fully removed
    assert False not in ['pgMag' not in i for i in fgs._present_calculated_mags]
    # check that gsc bands are all excluded from calculations
    assert fgs.survey == 'sdss'
    gsc_index = [i for i in fgs.band_dataframe['ABMag'].index if 'pgMag' in i]
    assert all(-999 == fgs.band_dataframe['ABMag'][gsc_index].values)

    # Test GSC only - set tmass and SDSS g, z, and i to missing (keep r).
    fgs = FGSCountrate(guide_star_id=gs_id, guider=guider)
    fgs.gsc_series = copy.copy(GSC_SERIES)
    fgs.gsc_series['tmassJMag'] = -999
    fgs.gsc_series['tmassHMag'] = -999
    fgs.gsc_series['tmassKsMag'] = -999
    fgs.gsc_series['SDSSgMag'] = -999
    fgs.gsc_series['SDSSzMag'] = -999
    fgs.gsc_series['SDSSiMag'] = -999
    fgs.calc_jhk_mag(fgs.gsc_series)
    _ = fgs.calc_fgs_cr_mag_and_err()
    # check sdss bands have been fully removed
    assert False not in ['SDSS' not in i for i in fgs._present_calculated_mags]
    # check that sdss bands are all excluded from calculations
    assert fgs.survey == 'gsc2'
    sdss_index = [i for i in fgs.band_dataframe['ABMag'].index if 'SDSS' in i]
    assert all(-999 == fgs.band_dataframe['ABMag'][sdss_index].values)

    # Test 2MASS only - set everything to missing except for tmass
    fgs = FGSCountrate(guide_star_id=gs_id, guider=guider)
    fgs.gsc_series = copy.copy(GSC_SERIES)
    fgs.gsc_series['SDSSgMag'] = -999
    fgs.gsc_series['SDSSzMag'] = -999
    fgs.gsc_series['SDSSiMag'] = -999
    fgs.gsc_series['SDSSrMag'] = -999
    fgs.gsc_series['SDSSuMag'] = -999
    fgs.gsc_series['JpgMag'] = -999
    fgs.gsc_series['FpgMag'] = -999
    fgs.gsc_series['NpgMag'] = -999
    fgs.calc_jhk_mag(fgs.gsc_series)
    _ = fgs.calc_fgs_cr_mag_and_err()
    # check sdss and gsc bands have been fully removed
    assert False not in ['SDSS' not in i for i in fgs._present_calculated_mags]
    assert False not in ['pgMag' not in i for i in fgs._present_calculated_mags]
    # check that sdss and gsc bands are all excluded from calculations
    assert fgs.survey == 'tmass'
    sdss_index = [i for i in fgs.band_dataframe['ABMag'].index if 'SDSS' in i]
    assert all(-999 == fgs.band_dataframe['ABMag'][sdss_index].values)
    gsc_index = [i for i in fgs.band_dataframe['ABMag'].index if 'pgMag' in i]
    assert all(-999 == fgs.band_dataframe['ABMag'][gsc_index].values)


def test_sdss_or_gsc_all_combinations():
    """ Full test of all combinations of bands to confirm that they correctly use only SDSS or only GSC
    values (or neither) when calculating the count rate and magnitude
    """
    gs_id = 'N13I000018'
    guider = 1

    # Iterate through every combination of present magnitudes
    for l in range(0, len(GSC_BAND_NAMES) + 1):
        for present_calculated_mags in itertools.combinations(GSC_BAND_NAMES, l):
            fgs = FGSCountrate(guide_star_id=gs_id, guider=guider)
            fgs.gsc_series = copy.copy(GSC_SERIES)

            # Set everything to -999 except for the present mags
            missing_mags = set(fgscountrate.fgs_countrate_core.GSC_BAND_NAMES) - set(present_calculated_mags)
            for ind in missing_mags:
                fgs.gsc_series.loc[ind] = -999

            # Run JHK and mag/cr calculation functions
            fgs.calc_jhk_mag(fgs.gsc_series)

            try:
                _ = fgs.calc_fgs_cr_mag_and_err()
            except ValueError as e:
                assert 'Cannot compute FGS countrate & magnitude for a guide star' in str(e)
                continue

            # Check survey
            if ('sdss' in fgs.j_convert_method) or ('sdss' in fgs.h_convert_method) or (
                    'sdss' in fgs.k_convert_method):
                assert fgs.survey == 'sdss'
            elif ('gsc2' in fgs.j_convert_method) or ('gsc2' in fgs.h_convert_method) or (
                    'gsc2' in fgs.k_convert_method):
                assert fgs.survey == 'gsc2'
            else:
                if len({'SDSSuMag', 'SDSSgMag', 'SDSSrMag', 'SDSSiMag', 'SDSSzMag'} & set(present_calculated_mags)) != 0:
                    assert fgs.survey == 'sdss', f'Present mags of {present_calculated_mags} not flagged as survey=sdss'
                elif len({'JpgMag', 'FpgMag', 'NpgMag'} & set(present_calculated_mags)) != 0:
                    assert fgs.survey == 'gsc2', f'Present mags of {present_calculated_mags} not flagged as survey=gsc2'
                elif {'tmassJMag', 'tmassHMag', 'tmassKsMag'} == set(present_calculated_mags):
                    assert fgs.survey == 'tmass', f'Present mags of {present_calculated_mags} not flagged as survey=tmass'

            # Check the band_dataframe indexes and values make sense
            if fgs.survey == 'sdss':
                # Check at least some SDSS values are included and no GSC2 values are included
                assert len({'SDSSuMag', 'SDSSgMag', 'SDSSrMag', 'SDSSiMag', 'SDSSzMag'} & set(fgs.band_dataframe.index.tolist())) != 0
                assert len({'JpgMag', 'FpgMag', 'NpgMag'} & set(fgs.band_dataframe.index.tolist())) == 0

                # Check the correct survey's bands not in present_calculated_mags are set to -999
                for band in ({'SDSSuMag', 'SDSSgMag', 'SDSSrMag', 'SDSSiMag', 'SDSSzMag'} & missing_mags):
                    assert fgs.band_dataframe['Signal'][band] == -999
                for band in ({'SDSSuMag', 'SDSSgMag', 'SDSSrMag', 'SDSSiMag', 'SDSSzMag'} & set(present_calculated_mags)):
                    assert fgs.band_dataframe['Signal'][band] != -999

            elif fgs.survey == 'gsc2':
                # Check at least some GSC2 values are included and no SDSS values are included
                assert len({'JpgMag', 'FpgMag', 'NpgMag'} & set(fgs.band_dataframe.index.tolist())) != 0
                assert len({'SDSSuMag', 'SDSSgMag', 'SDSSrMag', 'SDSSiMag', 'SDSSzMag'} & set(fgs.band_dataframe.index.tolist())) == 0

                # Check the correct survey's bands not in present_calculated_mags are set to -999
                for band in ({'JpgMag', 'FpgMag', 'NpgMag'} & missing_mags):
                    assert fgs.band_dataframe['Signal'][band] == -999
                for band in ({'JpgMag', 'FpgMag', 'NpgMag'} & set(present_calculated_mags)):
                    assert fgs.band_dataframe['Signal'][band] != -999

            elif fgs.survey == 'tmass':
                # Check there are only tmass values
                assert len({'tmassJMag', 'tmassHMag', 'tmassKsMag'} & set(fgs.band_dataframe.index.tolist())) != 0
                assert len({'SDSSuMag', 'SDSSgMag', 'SDSSrMag', 'SDSSiMag', 'SDSSzMag'} & set(fgs.band_dataframe.index.tolist())) == 0
                assert len({'JpgMag', 'FpgMag', 'NpgMag'} & set(fgs.band_dataframe.index.tolist())) == 0

                # Check that J and K must not be -999 (H can be)
                assert fgs.band_dataframe['Signal']['tmassJMag'] != -999
                assert fgs.band_dataframe['Signal']['tmassKsMag'] != -999


def test_errors():
    """Test errors are raised properly """

    gs_id = 'N13I000018'
    guider = 1

    # Test 1: data is missing too many bands - only has J and H (so cannot compute K)
    fgs = FGSCountrate(guide_star_id=gs_id, guider=guider)
    fgs.gsc_series = fgscountrate.utils.query_gsc(gs_id=gs_id, catalog='GSC242').iloc[0]

    fgs._present_calculated_mags = ['tmassJMag', 'tmassHMag']
    for index in set(fgscountrate.fgs_countrate_core.GSC_BAND_NAMES) - set(fgs._present_calculated_mags):
        fgs.gsc_series.loc[index] = -999
    fgs._all_calculated_mag_series = fgs.gsc_series.loc[fgscountrate.fgs_countrate_core.GSC_BAND_NAMES]
    fgs._all_calculated_mag_err_series = fgs.gsc_series.loc[[band+'Err' for band
                                                             in fgscountrate.fgs_countrate_core.GSC_BAND_NAMES]]
    fgs.survey = 'tmass'

    with pytest.raises(ValueError) as excinfo:
        fgs.calc_fgs_cr_mag_and_err()
    assert 'Cannot compute' in str(excinfo.value), 'Attempted to compute the FGS countrate & ' \
                                                   'magnitude despite only have J and H bands'

    # Test 2: Guider number is invalid
    guider = 3
    with pytest.raises(ValueError) as excinfo:
        fgs = FGSCountrate(guide_star_id=gs_id, guider=guider)
    assert '1 or 2' in str(excinfo.value), 'Allowed invalid guider number to pass'


def test_output_options():
    """
    Test the output options for calc_fgs_cr_mag_and_err()
    and _calc_fgs_cr_mag() are as expected
    """

    gs_id = 'N13I000018'
    guider = 2
    fgs = FGSCountrate(guide_star_id=gs_id, guider=guider)
    fgs.gsc_series = fgscountrate.utils.query_gsc(gs_id=gs_id, catalog='GSC242').iloc[0]
    fgs._present_calculated_mags = ['tmassJMag', 'tmassHMag', 'tmassKsMag', 'SDSSgMag', 'SDSSrMag', 'SDSSiMag']
    fgs._all_calculated_mag_series = fgs.gsc_series.loc[fgscountrate.GSC_BAND_NAMES]
    mag_err_list = [fgs.gsc_series[ind + 'Err'] for ind in fgs._all_calculated_mag_series.index]
    fgs._all_calculated_mag_err_series = pd.Series(mag_err_list, index=fgs._all_calculated_mag_series.index+'Err')
    fgs.survey = 'sdss'

    # Test output from calc_fgs_cr_mag_and_err()
    return_list = fgs.calc_fgs_cr_mag_and_err()
    assert len(return_list) == 4

    # Test output from _calc_fgs_cr_mag()
    band_series = fgs._all_calculated_mag_series
    guider_throughput = fgscountrate.fgs_countrate_core.THROUGHPUT_G2
    guider_gain = fgscountrate.fgs_countrate_core.CR_CONVERSION_G2

    # Case 1: Only Countrate
    return_list = fgs._calc_fgs_cr_mag(to_compute='countrate',
                                       band_series=band_series, guider_throughput=guider_throughput,
                                       guider_gain=guider_gain, return_dataframe=False)
    assert len(return_list) == 1
    assert np.isclose(return_list[0], fgs.fgs_countrate, 1e-5)

    # Case 2: Only Magnitude
    return_list = fgs._calc_fgs_cr_mag(to_compute='magnitude',
                                       band_series=band_series, guider_throughput=guider_throughput,
                                       guider_gain=guider_gain, return_dataframe=False)
    assert len(return_list) == 1
    assert np.isclose(return_list[0], fgs.fgs_magnitude, 1e-5)

    # Case 3: Both
    return_list = fgs._calc_fgs_cr_mag(to_compute='both',
                                       band_series=band_series, guider_throughput=guider_throughput,
                                       guider_gain=guider_gain, return_dataframe=False)
    assert len(return_list) == 2
    assert np.isclose(return_list[0], fgs.fgs_countrate, 1e-5)
    assert np.isclose(return_list[1], fgs.fgs_magnitude, 1e-5)

    # Case 4: Countrate + Dataframe Only
    return_list = fgs._calc_fgs_cr_mag(to_compute='countrate',
                                       band_series=band_series, guider_throughput=guider_throughput,
                                       guider_gain=guider_gain, return_dataframe=True)
    assert len(return_list) == 2
    assert np.isclose(return_list[0], fgs.fgs_countrate, 1e-5)
    np.testing.assert_array_almost_equal(return_list[1].values.flatten().tolist(),
                                         fgs.band_dataframe.values.flatten().tolist(), 5)

    # Case 5: Magnitude + Dataframe Only
    return_list = fgs._calc_fgs_cr_mag(to_compute='magnitude',
                                       band_series=band_series, guider_throughput=guider_throughput,
                                       guider_gain=guider_gain, return_dataframe=True)
    assert len(return_list) == 2
    assert np.isclose(return_list[0], fgs.fgs_magnitude, 1e-5)
    np.testing.assert_array_almost_equal(return_list[1].values.flatten().tolist(),
                                         fgs.band_dataframe.values.flatten().tolist(), 5)

    # Case 6: Countrate, Magnitude, and Dataframe
    return_list = fgs._calc_fgs_cr_mag(to_compute='both',
                                       band_series=band_series, guider_throughput=guider_throughput,
                                       guider_gain=guider_gain, return_dataframe=True)
    assert len(return_list) == 3
    assert np.isclose(return_list[0], fgs.fgs_countrate, 1e-5)
    assert np.isclose(return_list[1], fgs.fgs_magnitude, 1e-5)
    np.testing.assert_array_almost_equal(return_list[2].values.flatten().tolist(),
                                         fgs.band_dataframe.values.flatten().tolist(), 5)
