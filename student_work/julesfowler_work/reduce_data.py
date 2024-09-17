## -- IMPORTS

import glob

from astropy.io import fits
from astropy.stats import sigma_clip
import matplotlib.pyplot as plt
import numpy as np
from photutils import centroid_1dg, centroid_2dg


## -- FUNCTIONS

def calibrate_data(biases, darks, flats, science, writeout, plot):
    """ Main function to calibrate data from start to finish. 

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
    dark_image = make_dark(darks, bias_image)
    flat_image = make_flat(flats, dark_image)
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


def make_flat(flat_path, dark_image, writeout=None):
    """ Builds the dark subtracted flat image. 

    Parameters
    ----------
    flat_path: str
        Glob-able path to the flat fits files.
    dark_image : np.array
        Super diark in numpy array form.
    writeout : str or None
        Name of the fits file to write out, if at all.

    Returns
    -------
    flat : np.array
        Array of dark subtracted and time averaged flat data.
    """
   
    flat = median_combine(flat_path, norm_to_exposure=True, norm_to_1=True, subtract=dark_image)

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
        plt.clf()
        plt.imshow(reduced_science,
                vmin=np.median(reduced_science)-3*np.std(reduced_science),
                vmax=np.median(reduced_science)+3*np.std(reduced_science),
                cmap='bone')
        plt.savefig(plot)

    return reduced_science


def median_combine(path, norm_to_exposure=False, norm_to_1=False, subtract=None):
    """ Median combines a given set of images. Some boolean flags because flats
    are complicated.
    
    Parameters
    ----------
    path : str
        Glob-able path of files.
    norm_to_exposure : bool
        Whether or not to normalize by exposure time before median combining. 
    norm_to_1 : bool
        Whether or not to normalize to 1 before median combining.
    subtract : np.array or None
        Dark to subtract off before combining, or None.

    Returns
    -------
    median : np.array
        Array of median combined images.
    """

    files = glob.glob(path)
    size = np.shape(fits.getdata(files[0]))
    
    cube = np.zeros((size[0], size[1], len(files)))
    for index, fits_file in enumerate(files):

        with fits.open(fits_file) as hdu:
            data = hdu[0].data
            exp_time = hdu[0].header['EXPTIME']

        # Extra fiddling according to flags
        data = data/exp_time if norm_to_exposure else data
        data = data - subtract if subtract is not None else data
        sigma_max = np.max(sigma_clip(data, masked=False, axis=None))
        data = data/sigma_max if norm_to_1 else data
        
        cube[:, :, index] = data

    median = np.median(cube, axis=2)

    return median


def weighted_median_combine(path):
    """ Median combines a given set of images and weights them by their exposure time.
    This may be irrelevant now  but leaving it jic.

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

## -- pluto centroid night 1 falls at (496.90183365, 487.08766776)

## -- RUN
if __name__ == "__main__":
    
    # for pluto
    calibrate_data(biases='./pluto_data/biases/d*.fits', darks='./pluto_data/darks/night/d*.fits', 
                   flats='./pluto_data/sky_flats/d*fits', science='./pluto_data/science/d*.fits',
                   writeout='skyflats_pluto.fits', plot='skyflats_pluto.png')

    calibrate_data(biases='./pluto_data/biases/d*.fits', darks='./pluto_data/darks/night/d*.fits', 
                   flats='./pluto_data/dome_flats/d*fits', science='./pluto_data/science/d*.fits',
                   writeout='domeflats_pluto.fits', plot='domeflats_pluto.png')
    
    # for HR 
    # landolt B
    calibrate_data(biases='./hr_data/biases/d*.fits', darks='./hr_data/darks/d*.fits', 
                   flats='./hr_data/dome_flats/B_band/d*fits', science='./hr_data/science/landolt_B/d*.fits',
                   writeout='final_landolt_B.fits', plot='final_landolt_B.png')
    
    # landolt V
    calibrate_data(biases='./hr_data/biases/d*.fits', darks='./hr_data/darks/d*.fits', 
                   flats='./hr_data/dome_flats/V_band/d*fits', science='./hr_data/science/landolt_V/d*.fits',
                   writeout='final_landolt_V.fits', plot='final_landolt_V.png')
    
    # ngc B
    calibrate_data(biases='./hr_data/biases/d*.fits', darks='./hr_data/darks/d*.fits', 
                   flats='./hr_data/dome_flats/B_band/d*fits',
                   science='./hr_data/science/6819_B/d*.fits',
                   writeout='final_ngc6819_B.fits', plot='final_ngc6819_B.png')
    
    # ngc V
    calibrate_data(biases='./hr_data/biases/d*.fits', darks='./hr_data/darks/d*.fits', 
                   flats='./hr_data/dome_flats/V_band/d*fits',
                   science='./hr_data/science/6819_V/d*.fits',
                   writeout='final_ngc6819_V.fits', plot='final_ngc6819_V.png')
