import os
from pathlib import Path
import astropy.io.fits as fits
import numpy as np

read_path = Pathwrite_path = Path('SPICY_SE')
write_path = Pathwrite_path = Path('SPICY_SE_OUT')

files = [f for f in os.listdir(read_path) if f.endswith('.fits')]
files=sorted(files)
DETECT_THRESH=0.8
, 
for file in files:
    root = file.removesuffix('.fits')
    img = read_path / file
    chkimg = write_path / (root + '_SEG.fits')
    catlg = write_path / (root + '_CAT.fits')
    sexcmd1 = 'sextractor {} -c hosts.sex -CATALOG_NAME {} -CHECKIMAGE_NAME {} '.format(img, catlg, chkimg)
    sexcmd2 = '-CHECKIMAGE_TYPE SEGMENTATION -CATALOG_TYPE FITS_1.0 -DETECT_THRESH {}'.format(DETECT_THRESH)
    os.system(sexcmd1 + sexcmd2 + '>hosts_sex.log 2>&1')

# sextractor logfile in hosts_sex.log
# this script and provided sextractor config files - hosts.sex, hosts.param, default.conv - must be in the same folder

