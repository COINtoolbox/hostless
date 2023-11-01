import os
from pathlib import Path
import pandas as pd
import io
import gzip
import astropy.io.fits as fits
import numpy as np
from datetime import datetime
from datetime import timezone

read_path = Path('/media3/CRP7/hosts/data/SPICY/SPICY_CROSSMATCHED_2_ASEC_SMALL')
write_path = Path('SPICY_SE')

stamps = pd.read_parquet(read_path)

GAIN = 6.2       #e-/ADU
RDNOISE = 8      #e-
FULLWELL = 37000 #e-
SCALE = 1.0      #arcsec/pixel

field_name = ['SCI', 'TEMP', 'DIFF']

field_column = [
    'b:cutoutScience_stampData_small',
    'b:cutoutTemplate_stampData_small',
    'b:cutoutDifference_stampData_small'
]

for indx in range(len(stamps)):
    objid = stamps.iloc[indx]['objectId']
    times = stamps.iloc[indx]['i:jd']
    fwhm = stamps.iloc[indx]['i:fwhm']

    # these are to select only first and last stamp
    times_selected = [times[0],times[-1]]
    fwhm_selected = [fwhm[0],fwhm[-1]]
    stamp_selected = [0, len(times)-1]
    
    for field in range(len(field_name)):
        field_data = {}
        field_data[field] = stamps.iloc[indx][field_column[field]]        
        # for stamp_id in range(len(times)):  #if creating all stamps, not just selected,
        # use this instead of loop below and remove the "_selected" in all lines below.
        for stamp_id in [0,len(times_selected) -1]:
            stamp_time = times_selected[stamp_id]
            stamp_fwhm = fwhm_selected[stamp_id]
            stamp_num = stamp_selected[stamp_id]
            hdu_data_string = field_data[field][stamp_id]
            hdu = fits.open(gzip.open(io.BytesIO(hdu_data_string)))
            primary_hdu = hdu[0]
            hdr = primary_hdu.header
            hdr['OBJECT'] = (objid, 'string')
            hdr['IMTYPE'] = (field_name[field], 'string')
            hdr['DATE-OBS'] = (stamp_time, 'JD')
            hdr['STAMP'] = (str(stamp_num), 'integer')
            hdr['FWHM'] = (round(stamp_fwhm,2), 'pixel')
            hdr['SCALE'] = (SCALE, 'arcsec/pixel')
            hdr['GAIN'] = (GAIN, 'e-/ADU')
            hdr['RDNOISE'] = (RDNOISE, 'e-')
            hdr['FULLWELL'] = (FULLWELL, 'e-')
            created= datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%S%Z')
            hdr['HISTORY'] = ('Created ' + created)

            rootName  = objid + '_' + field_name[field] + '_' + str(stamp_num).zfill(4)
            StampFile = write_path  / str(rootName+'.fits')
            primary_hdu.writeto(StampFile, overwrite=True)
           

