"""
==================
FGS Countrate Core
==================

A module to find the expected FGS countrate for any JWST guide star

Code by Shannon Osborne. Contact at sosborne@stsci.edu
"""

from collections import OrderedDict
import io
import os
import requests
import sys

import numpy as np
import pandas as pd


class FGS_Countrate():
    """
    Class to support conversion from inputting a guide star ID and
    returning the expected FGS count rate. The main method in this
    class is get_fgs_countrate() which does that conversion.

    Parameters
    ----------
    guide_star_id : str
        The ID of the guide star of interest. This corresponds to the
        HST ID input in the Guide Star Catalog 2

    """

    def __init__(self, guide_star_id):
        self.id = guide_star_id

        self.data = None
        self.j_mag, self.h_mag, self.k_mag = None, None, None
        self.fgs_countrate = None

    def query_gsc(self):
        """
        Query the Guide Star Catalog 2.4.1 using the guide star ID.

        Returns
        -------
        data : pandas series
            A pd series containing the line from the GSC 2.4.1
            corresponding this the specific guide star ID

        """

        # Query GSC
        file_format = 'CSV'
        catalog = 'GSC241'
        url = 'http://gsss.stsci.edu/webservices/vo/CatalogSearch.aspx?' \
              'GSC2ID={}&FORMAT={}&CAT={}&'.format(self.id, file_format, catalog)
        request = requests.get(url).content

        # Read in data
        try:
            data_frame = pd.read_csv(io.StringIO(request.decode('utf-8')), skiprows=1)
            data_frame.replace(r'^\s+$', -999, regex=True, inplace=True)
        except pd.errors.EmptyDataError:
            raise NameError("This Guide Star ID does not exist in GSC2.4.1")

        # Check length of data table and turn it from a data frame to a series
        if len(data_frame) == 1:
            self.data = data_frame.iloc[0]
        else:
            # TODO: May do more fixing here later
            raise ValueError("This Guide Star ID points to multiple lines in GSC2.4.1")

        return self.data

    def get_fgs_countrate(self):
        """
        TBD
        """

        # Query GSC to get data
        data = self.query_gsc()

        # Convert to JHK magnitudes
        self.j_mag, self.h_mag, self.k_mag = self.convert_mag_to_jhk(data)

        # call convert_jhk_to_countrate
        return self.fgs_countrate

    def convert_mag_to_jhk(self, data):
        """
        Calculate the J, H, and K magnitudes of thee input data

        Parameters
        ----------
        data : pandas series
            A pd series containing the following stellar data:
            JpgMag, FpgMag, NpgMag, tmassJmag, tmassHmag,
            tmassKsMag, SDSSuMag, SDSSgMag, SDSSrMag,
            SDSSiMag, SDSSzMag

        Returns
        -------
        j, h, k : floats
            Returns the j, h, and k magnitudes for the star

        """

        # Pull the magnitudes section of the series
        l = ['JpgMag', 'FpgMag', 'NpgMag', 'tmassJmag', 'tmassHmag',
             'tmassKsMag', 'SDSSuMag', 'SDSSgMag', 'SDSSrMag',
             'SDSSiMag', 'SDSSzMag']
        mag_series = data.loc[l]

        # Rename variables for clarity
        mag_series = mag_series.rename({
            'JpgMag': 'B_Jmag',
            'FpgMag': 'R_Fmag',
            'NpgMag': 'I_Nmag',
            'tmassJmag': 'Jmag',
            'tmassHmag': 'Hmag',
            'tmassKsMag': 'Kmag',
            'SDSSuMag': 'umag',
            'SDSSgMag': 'gmag',
            'SDSSrMag': 'rmag',
            'SDSSiMag': 'imag',
            'SDSSzMag': 'zmag'})

        # List of the magnitude names that are not fill values in the series
        present_mags = list(mag_series[mag_series != -999].index)

        # Create a list of magnitudes that are important in determining conversion method
        present_mags_short = [e for e in present_mags if e in ['B_Jmag', 'R_Fmag', 'I_Nmag', 'Jmag',
                                                               'gmag', 'imag', 'zmag']]

        # Dictionary of convert methods
        switcher = OrderedDict([
            ('Jmag', 'convert_1'),
            ('gmag, zmag', 'convert_2'),
            ('gmag', 'convert_3'),
            ('imag, zmag', 'convert_4'),
            ('imag', 'convert_5'),
            ('B_Jmag, I_Nmag', 'convert_6'),
            ('B_Jmag, R_Fmag', 'convert_7'),
            ('B_Jmag', 'convert_8'),
            ('R_Fmag, I_Nmag', 'convert_9'),
            ('R_Fmag', 'convert_10'),
            ('', 'convert_11noinfo'),
        ])

        # Pull the first entry in the OrderedDict that matches what values are present.
        for key, value in switcher.items():
            key_list = key.split(', ')
            if set(key_list).issubset(present_mags_short):
                name = key
                break
            else:
                name = ''  # there isn't enough information to follow the flowchart

        # Get the method
        method_name = switcher.get(name, "nothing")
        method = getattr(self, method_name, lambda: "Invalid")

        return method()

    # These may go in another file for cleanliness. We'll see how long they are
    def convert_1(self):
        """No conversion needed. Input is already in j,h,k bands"""
        return j, h, k

    def convert_2(self):
        return j, h, k

    def convert_3(self):
        return j, h, k

    def convert_11noinfo(self):
        raise ValueError('There is not enough information on this guide star to get its J, H and K data')


    def convert_jhk_to_countrate(self):
        # ...
        # self.fgs_countrate = value
        return self.fgs_countrate