import itertools

import pytest

from fgscountrate.fgs_countrate_core import FGSCountrate
from fgscountrate import utils


def test_successful_query():
    """Test a query that should be successful"""
    gs_id = 'N13I000018'
    dataframe = utils.query_gsc(gs_id=gs_id, catalog='GSC242')

    assert len(dataframe) == 1

    gs_ra = dataframe['ra'][0]
    gs_dec = dataframe['dec'][0]

    assert pytest.approx(gs_ra, 1e-5) == 273.206729760604
    assert pytest.approx(gs_dec, 1e-5) == 65.5335149359777


# TODO: Get an example ID for this test
@pytest.mark.skip(reason="Don't have an example ID for this test yet")
def test_multiple_line_output():
    """Check that an error occurs when an ID corresponds to multiple lines in the GSC"""
    gs_id = 'N94D006564'

    with pytest.raises(ValueError) as excinfo:
        utils.query_gsc(gs_id=gs_id)
    assert 'multiple lines' in str(excinfo.value), 'GSC only finding one line for an ID which has multiple lines'


def test_fake_guidestar_id():
    """Check errors built into querying method"""

    # Check an input guide star ID that isn't found
    with pytest.raises(NameError) as excinfo:
        gs_id = 'fakeid'
        utils.query_gsc(gs_id=gs_id, catalog='GSC242')
    assert 'No guide stars match' in str(excinfo.value), 'Fake Guide Star ID incorrectly found in catalog'


def test_coordinate_combinations():
    """
    Test different possible combinations of coordinates to query the GSC

    These inputs were checked for GSC2.4.1 and should return guide stars. A
    change to the guide star catalog could cause this test to throw an error
    if the guide star entries have changed
    """

    input_dict = {
        'gs_id': 'N13I000018',
        'ra': 273.206729584718,
        'dec': 65.5335161247318,
        'cone_radius': 0.1,
        'minra': 273.2,
        'maxra': 273.21,
        'mindec': 65.5,
        'maxdec': 65.55,

    }

    for L in range(0, len(input_dict) + 1):
        for subset in itertools.combinations(list(input_dict.keys()), L):

            # Check only 1 coordinate specification/a full set of coordinates is being used
            if set(subset) not in [set(['gs_id']), set(['ra', 'dec']), set(['ra', 'dec', 'cone_radius']),
                                   set(['minra', 'maxra', 'mindec', 'maxdec'])]:
                with pytest.raises(ValueError):
                    subset_dict = dict((k, input_dict[k]) for k in subset)
                    utils.query_gsc(**subset_dict, catalog='GSC242')

            # Check that the remaining combinations should successfully run
            else:
                print(subset)
                subset_dict = dict((k, input_dict[k]) for k in subset)
                utils.query_gsc(**subset_dict, catalog='GSC242')
