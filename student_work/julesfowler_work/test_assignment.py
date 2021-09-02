""" Quick script to fufill the first test python assignment. """

## -- IMPORTS

import matplotlib.pyplot as plt
import numpy as np

from astropy.io import fits
from photutils import centroid_1dg, centroid_2dg

## -- FUNCTIONS

def crop_around_max(image, region_size, name_out='test2'):
    """ Crops a region_size x region_size box around the brightest
    point in the image. (I did some interactive testing to make
    sure that brightest point wasn't an edge artifact before 
    letting a simple max do the trick.)

    Parameters
    ----------
    image : np.array 
        Array with image data.
    region_size : int
        Size of the square to scoop out.
    name_out : str, optional
        Name of the file to write out with the cropped region.

    Returns
    -------
    new_region : np.array
        Array of cropped image data.
    """
    
    # Find the brightest point and rip off the trappings
    max_value = np.nanmax(image)
    m, n = np.where(image == max_value)
    m, n = m[0], n[0]
    
    new_region = image[m-int(region_size/2):m+int(region_size/2), n-int(region_size/2):n+int(region_size/2)]
    
    if name_out is not None:
        hdu = fits.PrimaryHDU(new_region)
        hdu_list = fits.HDUList([hdu])
        hdu_list.writeto(f'{name_out}.fits', overwrite=True)
    return new_region


def find_centroid(cropped_image, method='compare', plot=True):
    """ Finds the centroid in an image. 

    Parameters
    ----------
    cropped_image : np.array
        Array of image data.
    method : str, optional
        Photutils method to take the centroid.
    plot : bool
        Whether or not to plot centroid.
    """ 

    if method == '1dg':
        photutils_method = [(centroid_1dg, method)]
    elif method == '2dg':
        photutils_method = [(centroid_2dg, method)]
    elif method == 'compare':
        photutils_method = [(centroid_1dg, '1dg'), (centroid_2dg, '2dg')]
    else:
        raise NotImplementedError("Sorry, I can't run the centroid method you want. Try '1dg', '2dg', or 'compare'.")

    plt.imshow(cropped_image, cmap='bone')
    colors=['cyan', 'teal']
    for index, photometry in enumerate(photutils_method):
        m, n = photometry[0](cropped_image)
        print(f'Found a centroid position of ({m},{n}) with the {photometry[1]} method.')
        plt.scatter(m, n, facecolors='none', edgecolors=colors[index], s=500)

    if plot:
        plt.savefig('centroid.png')


def main():
    """ Main function to scoop up file, crop it, and run photometry."""

    image_file = '../../python/test_assignment/test.fits'
    try:
        data = fits.getdata(image_file)
    except FileNotFoundError:
        raise FileNotFoundError("No test.fits exists in the local directory. Are we in the right place?")
    
    cropped_image = crop_around_max(data, region_size=40)
    find_centroid(cropped_image)


## -- RUN

if __name__ == "__main__":
    main()




