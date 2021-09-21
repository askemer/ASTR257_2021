## -- IMPORTS

import glob

from astropy.io import fits
from astropy.stats import sigma_clip
import matplotlib.pyplot as plt
import numpy as np
from photutils import centroid_1dg, centroid_2dg


## -- FUNCTIONS

def calibrate_data(biases='./biases/d*.fits', darks='./darks/night/d*.fits', 
                   flats='./dome_flats/d*fits', science='./science/d*.fits',
                   writeout='final_pluto.fits', plot='final_pluto.png'):
    """ Main function to calibrate data from start to finish. 
    Note that with the default paths this is meant to run in a particuarly
    data directory structure. Future iterations will probably change that.

    Parameters
    ----------
    biases: str
        Glob-able path to bias fits files.
    darks: str
        Glob-able path to dark fits files.
    flats: str
        Glob-able path to flat fits files. (Will work with sky or dome flats.
    science: str
        Glob-able path to science fits files.
    writeout: str
        Name to write out the final calibrated pluto science image to.
    plot: str
        Name to write out the final calibrated pluto plot to.
    
    Returns
    -------
    science_image : np.array
        Array of reduced science data.
    """

    bias_image = make_bias(biases)
    dark_image = make_dark(darks, bias_image, writeout='dark.fits')
    flat_image = make_flat(flats, dark_image, writeout='flat.fits')
    science_image = make_science_image(science, bias_image, dark_image, 
            flat_image, writeout=writeout, plot=plot)

    return science_image


def make_bias(bias_path):
    """ Quick one-liner to return the super bias in a modular way."""
    return median_combine(bias_path)

def make_dark(dark_path, bias_image, writeout=None):
    """ Builds the bias subtracted dark image.
    
    Parameters
    ----------
    dark_path : str
        Glob-able path to dark fits files.
    bias_image : np.array
        Super bias in numpy array form.
    writeout : str or None
        Name of the fits file to write out, if at all.
    
    Returns
    -------
    super_dark_per_second : np.array
        Bias subtracted + time averaged super dark.
    """

    super_dark = median_combine(dark_path)
    super_dark_subtracted = super_dark - bias_image
    
    dark_exp_time = fits.getheader(glob.glob(dark_path)[0])['EXPTIME']
    super_dark_per_second = super_dark_subtracted/dark_exp_time 
    
    if writeout is not None:
        hdu_list = fits.HDUList([fits.PrimaryHDU(super_dark_per_second)])
        hdu_list.writeto(writeout, overwrite=True)

    return super_dark_per_second


def make_flat(flat_path, dark_image, norm_to_1=True, writeout=None):
    """ Builds the dark subtracted flat image. 

    Parameters
    ----------
    flat_path: str
        Glob-able path to the flat fits files.
    dark_image : np.array
        Super diark in numpy array form.
    norm_to_1: bool
        Whether or not to normalize the flat to 1.
    writeout : str or None
        Name of the fits file to write out, if at all.

    Returns
    -------
    flat : np.array
        Array of dark subtracted and time averaged flat data.
    """
    super_flat = weighted_median_combine(flat_path)
    super_flat_subtracted = super_flat - dark_image

    if norm_to_1:
        sigma_max = np.max(sigma_clip(super_flat_subtracted, masked=False, axis=None))
        flat = super_flat_subtracted/sigma_max
    else:
        flat = super_flat_subtracted

    if writeout is not None:
        hdu_list = fits.HDUList([fits.PrimaryHDU(flat)])
        hdu_list.writeto(writeout, overwrite=True)

    return flat


def make_science_image(science_image_path, bias_image, dark_image, flat_image, writeout=None, plot=None):
    """ Builds the final reduced science image. 

    Parameters
    ----------
    science_image_path: str
        Glob-able path to the science fits files.
    bias_image : np.array
        Array of super bias data.
    dark_image : np.array
        Array of super dark data.
    flat_image : np.array
        Dark subtracted and time averaged flat data.
    writeout : str or None
        Name of the fits file to write out, if at all.
    plot : str or None
        Name of the plot to write out, if at all.

    Returns
    -------
    reduced_science : np.array
        Array of final reduced science image.
    """

    super_science = median_combine(science_image_path)
    super_science_subtracted = super_science - bias_image
    
    science_exp_time = fits.getheader(glob.glob(science_image_path)[0])['EXPTIME']
    super_science_per_second = super_science_subtracted/science_exp_time

    reduced_science = (super_science_per_second - dark_image)/flat_image

    if writeout is not None:
        hdu_list = fits.HDUList([fits.PrimaryHDU(reduced_science)])
        hdu_list.writeto(writeout, overwrite=True)

    if plot is not None:
        plt.imshow(reduced_science,
                vmin=np.median(reduced_science)-3*np.std(reduced_science),
                vmax=np.median(reduced_science)+3*np.std(reduced_science),
                cmap='bone')
        plt.savefig(plot)

    return reduced_science


def median_combine(path):
    """ Median combines a given set of images.
    
    Parameters
    ----------
    path : str
        Glob-able path of files.

    Returns
    -------
    median : np.array
        Array of median combined images.
    """

    files = glob.glob(path)
    size = np.shape(fits.getdata(files[0]))
    
    cube = np.zeros((size[0], size[1], len(files)))
    for index, fits_file in enumerate(files):
        cube[:, :, index] = fits.getdata(fits_file)

    median = np.median(cube, axis=2)

    return median


def weighted_median_combine(path):
    """ Median combines a given set of images and weights them by their exposure time.

    Parameters
    ----------
    path : str
        Glob-able path of files.

    Returns
    -------
    median : np.array
        Array of median combined + weighted images.
    """

    files = glob.glob(path)
    size = np.shape(fits.getdata(files[0]))

    cube = np.zeros((size[0], size[1], len(files)))
    for index, fits_file in enumerate(files):
        with fits.open(fits_file) as hdu:
            data = hdu[0].data
            exp_time = hdu[0].header['EXPTIME']
            cube[:, :, index] = data/exp_time if exp_time != 0 else data

    median = np.median(cube, axis=2)
    return median

## -- RUN
if __name__ == "__main__":
    calibrate_data()
