import io

import numpy as np
import pandas as pd
import requests


def convert_to_abmag(value, name):
    """
    Convert magnitude to AB magnitude

    Parameters
    ----------
    value : float
        Value of the band
    name : str
        Name of the band as stated in the GSC column name.
        Options are: 2MASS: tmassJmag, tmassHmag, tmassKsMag
        SDSS: SDSSgMag, SDSSiMag, SDSSzMag
        GSC: JpgMag, FpgMag, IpgMag

    """

    mag_constants = {
        'tmassJmag': 0.90,
        'tmassHmag': 1.37,
        'tmassKsMag': 1.85,
        'SDSSuMag': 0.0,
        'SDSSgMag': 0.0,
        'SDSSrMag': 0.0,
        'SDSSiMag': 0.0,
        'SDSSzMag': 0.0,
        'JpgMag': -0.055,
        'FpgMag': 0.24,
        'NpgMag': 0.48,
    }

    abmag = value + mag_constants[name]

    return abmag
#GSC2_and_2MASS_and_SDSS_RA273.19deg_Dec65.54deg

def query_gsc(gs_id=None, ra=None, dec=None, cone_radius=None, minra=None, maxra=None,
              mindec=None, maxdec=None, catalog=None):
    """
    Query the Guide Star Catalog using one of 4 query options:
        1) Guide Star ID
        2) Exact RA & DEC
        3) Centering RA & DEC and a cone search radius
        4) Min/Max RA and Min/Max DEC for a bounding box search
    Only pass in values for one of the 4 query options and pass all
    required values

    Parameters
    ----------
    gs_id : str
        The ID of the guide star of interest. This corresponds to the
        HST ID input in the Guide Star Catalog 2
    ra : float
        The right ascension in degrees of the target or catalog sub-section
        to be retrieved.
    dec : float
        The declination in degrees of the target or catalog sub-section
        to be retrieved.
    cone_radius : float
        Cone search radius in degrees.
    minra : float
        Minimum right ascension in degrees for box search.
    maxra : float
        Maximum right ascension in degrees for box search.
    mindec : float
        Minimum declination in degrees for box search.
    maxdec : float
        Maximum declination in degrees for box search.
    catalog : str
        There are 5 different GSC2 versions available. Default is GSC241
        Call GSC23 to access GSC2.3.4
        Call GSC240 to access GSC2.4.0
        Call GSC241 to access GSC2.4.1.1
        Call GSC2412 to access GSC2.4.1.2
        Call GSC2420 to access GSC2.4.2

    Returns
    -------
    data_frame : Pandas DataFrame
        A pd dataframe containing the row(s) from the specified catalog
        corresponding to the requested GS ID, coordinates, &/or area

    """

    # Set file format and default catalog
    file_format = 'CSV'
    if catalog is None:
        catalog = 'GSC241'

    # Check only 1 coordinate specification is being used AND the coordinate specification chosen is complete
    method_list = [any([gs_id]), any([ra, dec, cone_radius]), any([minra, maxra, mindec, maxdec])]
    complete_list = [all([gs_id]), all([ra, dec]) or all([ra, dec, cone_radius]), all([minra, maxra, mindec, maxdec])]
    if method_list.count(True) != 1:
        raise ValueError("You may only specify coordinates using one method.")
    if complete_list.count(True) != 1:
        raise ValueError("You must specify a full set of coordinates for your chosen method.")

    # Set URL
    url = 'http://gsss.stsci.edu/webservices/vo/CatalogSearch.aspx?'
    if gs_id is not None:
        url = url + 'GSC2ID={}&'.format(gs_id)
    if ra is not None:
        url = url + 'RA={}&'.format(ra)
    if dec is not None:
        url = url + 'DEC={}&'.format(dec)
    if cone_radius is not None:
        url = url + 'SR={}&'.format(cone_radius)
    if minra is not None:
        url = url + 'BBOX={}%2c'.format(minra)
    if mindec is not None:
        url = url + '{}%2c'.format(mindec)
    if maxra is not None:
        url = url + '{}%2c'.format(maxra)
    if maxdec is not None:
        url = url + '{}&'.format(maxdec)
    url = url + 'FORMAT={}&CAT={}'.format(file_format, catalog)

    # Query data
    request = requests.get(url).content

    # Read data into pandas
    try:
        data_frame = pd.read_csv(io.StringIO(request.decode('utf-8')), skiprows=1, na_values=[' '])
        data_frame.replace(np.nan, -999, regex=True, inplace=True)
    except pd.errors.EmptyDataError:
        raise NameError("No guide stars match these requirements in catalog {}".format(catalog))

    return data_frame
