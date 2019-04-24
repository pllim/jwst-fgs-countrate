JWST FGS countrate estimation
-----------------------------

[![STScI](https://img.shields.io/badge/powered%20by-STScI-blue.svg?colorA=707170&colorB=3e8ddd&style=flat)](http://www.stsci.edu)


jwst-fgs-countrate converts guide star information from different catalogs into the expected countrate for the JWST Fine Guidance Sensor (FGS). This is done by querying the Guide Star Catalog (containing GSC2, 2MASS, SDSS, and eventually GAIA data) based on the input of the guide star’s HST ID number and using information from that catalog to calculate a countrate.

Collaborators include Shannon Osborne ([@shanosborne](https://github.com/shanosborne)), Sherie Holfeltz ([@stholfeltz](https://github.com/stholfeltz)), and Pierre Chayer.


Example Use
-----------

To compute the FGS countrate given a guide star ID:
```
fgs = FGS_Countrate(guide_star_id='N13I000018')
fgs_countrate = get_fgs_countrate()
```

To get the guide star data from the GSC2.4.1:
```
# After doing the countrate calculation
fgs = FGS_Countrate(guide_star_id='N13I000018')
fgs_countrate = get_fgs_countrate()
fgs.data

# Or independent of the countrate calculation
fgs = FGS_Countrate(guide_star_id='N13I000018')
data = fgs.query_gsc()
```

To learn which conversions were used when calculating the J, H, and K magnitudes
```
fgs.j_convert_method
fgs.h_convert_method
fgs.k_convert_method
```

License
-------

This project is Copyright (c) Space Telescope Science Institute and licensed under
the terms of the Aura license. This package is based upon
the [Astropy package template](https://github.com/astropy/package-template)
which is licensed under the BSD 3-clause licence. See the licenses folder for
more information.
