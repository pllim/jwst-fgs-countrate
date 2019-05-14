import pytest

from fgscountrate.fgs_countrate_core import FGSCountrate
from fgscountrate import utils


def test_successful_query():
    """Test a query that should be successful"""
    id = 'N13I000018'
    fgs = FGSCountrate(guide_star_id=id, guider=1)
    dataframe = utils.query_gsc(gs_id=fgs.id, catalog='GSC241')

    assert len(dataframe) == 1

    gs_ra = dataframe['ra'][0]
    gs_dec = dataframe['dec'][0]

    assert pytest.approx(gs_ra, 1e-5) == 273.206729760604
    assert pytest.approx(gs_dec, 1e-5) == 65.5335149359777


def test_fake_id():
    """Check errors built into querying method"""
    id = 'fake_id'

    with pytest.raises(NameError) as excinfo:
        dataframe = utils.query_gsc(gs_id=id, catalog='GSC241')
    assert 'No guide stars match' in str(excinfo.value), 'Fake Guide Star ID incorrectly found in catalog'


@pytest.mark.skip(reason="Don't have an example ID for this test yet")
def test_multiple_line_output():
    """Check that an error occurs when an ID corresponds to multiple lines in the GSC"""
    id = 'TBD'
    fgs = FGSCountrate(guide_star_id=id, guider=2)

    with pytest.raises(ValueError) as excinfo:
        fgs.get_fgs_countrate()
    assert 'multiple lines' in str(excinfo.value), 'GSC only finding one line for an ID which has multiple lines'
