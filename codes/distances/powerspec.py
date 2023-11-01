import os
import glob
import numpy as np
import pandas as pd
from astropy.io import fits
from scipy.stats import binned_statistic
import argparse

parser = argparse.ArgumentParser(description="Fits power spectrum to stamps (originally of size 63x63) in 3 size bins (7x7, 15x15, 29x29).")
parser.add_argument("--candidates", nargs='+', help="Name or list of names of hostless candidates.", required=True)
args = parser.parse_args()

N = 1000
full = 63
size_dict = {'small': 7, 'medium': 15, 'large': 29}

small_len = len(np.arange(0.5, size_dict['small']//2+1, 1.))-1
med_len = len(np.arange(0.5, size_dict['medium']//2+1, 1.))-1
large_len = len(np.arange(0.5, size_dict['large']//2+1, 1.))-1

for objectID in args.candidates:
    print(f"Starting {objectID}")

    os.makedirs(f"/media3/hayes/crp7_hostless/pspec/{objectID}")
    real_dict = {'SCI': {'small': np.zeros(small_len), 'medium': np.zeros(med_len), 'large': np.zeros(large_len)},
                 'TEMP': {'small': np.zeros(small_len), 'medium': np.zeros(med_len), 'large': np.zeros(large_len)}}
    output_dict = {'SCI': {'small': np.zeros((N, small_len)), 'medium': np.zeros((N, med_len)), 'large': np.zeros((N, large_len))},
                   'TEMP': {'small': np.zeros((N, small_len)), 'medium': np.zeros((N, med_len)), 'large': np.zeros((N, large_len))}}

    for type in ['SCI', 'TEMP']:

        seg_files = glob.glob(f"/media3/CRP7/hosts/data/EXTRAGALACTIC/sexoutput/{objectID}_{type}_????_SEG.fits")
        im_files = glob.glob(f"/media3/CRP7/hosts/data/EXTRAGALACTIC/sexoutput/{objectID}_{type}_????.fits")
        seg = fits.open(seg_files[-1])
        image = fits.open(im_files[-1])

        mask = np.logical_or(seg[0].data > 0, np.where(np.isnan(image[0].data), True, False))
        for_filling = np.random.normal(np.median(image[0].data[~mask]), np.std(image[0].data[~mask]), (63, 63))
        for_filling = np.where(mask, for_filling, 0)
        to_fill = np.where(mask, 0, image[0].data)
        data = to_fill + for_filling
        
        for n in range(N):
            copy = np.copy(data)
            copy = copy.reshape(63*63)
            np.random.shuffle(copy)
            copy = copy.reshape((63, 63))

            for s, size in enumerate(size_dict.keys()):
                start = int((full - size_dict[size])/2)
                stop = int((full + size_dict[size])/2)

                if n == 0:
                    data_resized = data[start : stop, start : stop]
                    fourier_image = np.fft.fftn(data_resized)
                    fourier_amplitudes = np.abs(fourier_image)**2
                    kfreq = np.fft.fftfreq(size_dict[size]) * size_dict[size]
                    kfreq2D = np.meshgrid(kfreq, kfreq)
                    knrm = np.sqrt(kfreq2D[0]**2 + kfreq2D[1]**2)
                    knrm = knrm.flatten()
                    fourier_amplitudes = fourier_amplitudes.flatten()
                    kbins = np.arange(0.5, size_dict[size]//2+1, 1.)
                    kvals = 0.5 * (kbins[1:] + kbins[:-1])
                    Abins, _, _ = binned_statistic(knrm, fourier_amplitudes, statistic = "mean", bins = kbins)
                    Abins *= np.pi * (kbins[1:]**2 - kbins[:-1]**2)

                    real_dict[type][size] = Abins

                copy_resized = copy[start : stop, start : stop]

                fourier_image = np.fft.fftn(copy_resized)
                fourier_amplitudes = np.abs(fourier_image)**2
                kfreq = np.fft.fftfreq(size_dict[size]) * size_dict[size]
                kfreq2D = np.meshgrid(kfreq, kfreq)
                knrm = np.sqrt(kfreq2D[0]**2 + kfreq2D[1]**2)
                knrm = knrm.flatten()
                fourier_amplitudes = fourier_amplitudes.flatten()
                kbins = np.arange(0.5, size_dict[size]//2+1, 1.)
                kvals = 0.5 * (kbins[1:] + kbins[:-1])
                Abins, _, _ = binned_statistic(knrm, fourier_amplitudes, statistic = "mean", bins = kbins)
                Abins *= np.pi * (kbins[1:]**2 - kbins[:-1]**2)

                output_dict[type][size][n] = Abins

    for type in ['SCI', 'TEMP']:
        for size in size_dict.keys():
            array = pd.DataFrame(output_dict[type][size])
            array.to_csv(f"/media3/hayes/crp7_hostless/pspec/{objectID}/{type}_{size}_shuffled.csv")

            array = pd.DataFrame(real_dict[type][size])
            array.to_csv(f"/media3/hayes/crp7_hostless/pspec/{objectID}/{type}_{size}.csv")
        
    print(f"Done {objectID}")

