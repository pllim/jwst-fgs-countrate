import pytest

from fgscountrate.fgs_countrate_core import FGS_Countrate


def test_successful_query():
    """Test a query that should be successful"""
    id = 'N13I000018'
    fgs = FGS_Countrate(guide_star_id=id)
    data_series = fgs.query_gsc()

    gs_ra = data_series.loc['ra']
    gs_dec = data_series.loc['dec']

    assert pytest.approx(gs_ra, 1e-5) == 273.206729760604
    assert pytest.approx(gs_dec, 1e-5) == 65.5335149359777


def test_fake_id():
    """Check errors built into querying method"""
    id = 'fake_id'
    fgs = FGS_Countrate(guide_star_id=id)

    with pytest.raises(NameError) as excinfo:
        data_series = fgs.query_gsc()
    assert 'does not exist' in str(excinfo.value), 'Fake Guide Star ID incorrectly found in catalog'


@pytest.mark.skip(reason="Don't have an example ID for this test yet")
def test_multiple_line_output():
    """Check that an error occurs when an ID corresponds to multiple lines in the GSC"""
    id = 'TBD'
    fgs = FGS_Countrate(guide_star_id=id)

    with pytest.raises(ValueError) as excinfo:
        data_series = fgs.query_gsc()
    assert 'multiple lines' in str(excinfo.value), 'GSC only finding one line for an ID which has multiple lines'



