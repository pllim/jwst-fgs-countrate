"""
==================
FGS Countrate Core
==================

A module to find the expected FGS countrate for any JWST guide star

Code by Shannon Osborne. Contact at sosborne@stsci.edu
"""

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
        data : pandas data table
            A length = 1 table containing the line from the GSC 2.4.1
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
            self.data = pd.read_csv(io.StringIO(request.decode('utf-8')), skiprows=1)
            self.data.replace(r'^\s+$', -999, regex=True, inplace=True)
        except pd.errors.EmptyDataError:
            raise NameError("This Guide Star ID does not exist in GSC2.4.1")

        # Check length of data table
        if len(self.data) != 1:
            # TODO: May do more fixing here later
            raise ValueError("This Guide Star ID points to multiple lines in GSC2.4.1")

        # TODO: Parse data? TBD

        return self.data

    def get_fgs_countrate(self):
        # call convert_id_to_jhk
        # call convert_jhk_to_countrate
        return self.fgs_countrate

    def convert_id_to_jhk(self):
        # call query_gsc
        # ...
        # self.j_mag, self.h_mag, self.k_mag = value
        return self.j_mag, self.h_mag, self.k_mag

    # These may go in another file for cleanliness. We'll see how long they are
    def _case1(self):
        """No conversion needed. Input is already in j,h,k bands"""
        return j, h, k
    def _case2(self):
        return j, h, k
    def _case3(self):
        return j, h, k
    def _case4(self):
        return j, h, k


    def convert_jhk_to_countrate(self):
        # ...
        # self.fgs_countrate = value
        return self.fgs_countrate