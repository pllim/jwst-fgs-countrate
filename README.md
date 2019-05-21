JWST FGS countrate estimation
-----------------------------

[![STScI](https://img.shields.io/badge/powered%20by-STScI-blue.svg?colorA=707170&colorB=3e8ddd&style=flat)](http://www.stsci.edu)


jwst-fgs-countrate converts guide star information from different catalogs into the expected countrate for the JWST Fine Guidance Sensor (FGS). This is done by querying the Guide Star Catalog (containing GSC2, 2MASS, SDSS, and eventually GAIA data) based on the input of the guide star’s HST ID number and using information from that catalog to calculate a countrate.

Collaborators include Shannon Osborne ([@shanosborne](https://github.com/shanosborne)), Sherie Holfeltz ([@stholfeltz](https://github.com/stholfeltz)), and Pierre Chayer.


Example Use
-----------

To compute the FGS countrate and magnitude given a guide star ID:
```
import fgscountrate
fgs = fgscountrate.FGS_Countrate(guide_star_id='N13I000018', guider=1)
cr, cr_err, mag, mag_err = fgs.query_fgs_countrate_magnitude()
```

The countrate and magnitude calculations are also saved in attributes
```
# Countrate and Magnitude attributes
fgs.countrate
fgs.countrate_err
fgs.magnitude
fgs.magnitude_err
```

To get the guide star data from the guide star catalog
```
# After doing the countrate calculation, access the gsc_series attribute
import fgscountrate
fgs = FGS_Countrate(guide_star_id='N13I000018', guider=1)
fgs_countrate = fgs.query_fgs_countrate_magnitude()
fgs.gsc_series

# Or independent of the countrate calculation, use query_gsc()
from fgscountrate.utils import query_gsc
data = query_gsc(gs_id='N13I000018')
```

To learn which conversions were used when calculating the J, H, and K magnitudes
```
# Conversion method attributes
fgs.j_convert_method
fgs.h_convert_method
fgs.k_convert_method
```

Or to access any of the individual JHK conversion equations
```
# Using already known data
from fgscountrate import conversions
g_mag = 14.5
z_mag = 15.2
j_mag = conversions.convert_sdssgz_to_jhk(data=(g_mag, z_mag), output_mag='J')

# Using data queried from the GSC using query_gsc()
from fgscountrate.utils import query_gsc
dataframe = query_gsc(gs_id='N13I000018')
dataseries = dataframe.iloc[0]
j_mag = conversions.convert_sdssgz_to_jhk(data=dataseries, output_mag='J')
```

To get the data from the step by step calculations for the FGS countrate and magnitude
```
# Returns a dataframe
fgs.band_dataframe
```

License
-------

This project is Copyright (c) Space Telescope Science Institute and licensed under
the terms of the Aura license. This package is based upon
the [Astropy package template](https://github.com/astropy/package-template)
which is licensed under the BSD 3-clause licence. See the licenses folder for
more information.
