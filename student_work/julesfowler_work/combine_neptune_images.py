### -- NEPTUNE AO ANALYSIS
"""Holds the code for the Neptune AO lab. """

## -- IMPORTS
from astropy.io import fits
from image_registration import chi2_shift
from image_registration.fft_tools import shift
import numpy as np

## -- FUNCTIONS

def build_super_images(files, register=True):

    data = np.array([fits.getdata(img_file) for img_file in files])
    median_img = np.median(data, axis=0)

    base_img = data[0] - median_img
    corrected_images =  np.zeros_like(data)
    corrected_images[0] = base_img
    for i in range(1, len(data)):
        if register:
            corrected_image = register_image(base_img, data[i]-median_img)
        else:
            corrected_image = data[i] - median_img
        corrected_images[i] = corrected_image

    super_image = np.median(corrected_images, axis=0)
    return super_image, corrected_images

def register_image(original, image_to_correct):

    xoff, yoff, exoff, eyoff = chi2_shift(original, image_to_correct, np.std(original))
    corrected_image = shift.shiftnd(image_to_correct, (-yoff, -xoff))

    return corrected_image

## -- GO

if __name__ == "__main__":

    # H 
    h_files = glob.glob('h/s*.fits')
    super_image_h, _ = build_super_images(h_files)
    
    # J
    j_files = glob.glob('j/s*.fits')
    super_image_j, _ = build_super_images(j_files)

    # CH4
    # this one really struggled with image registration
    # so I'm just picking one good one and applying a median subtraction
    ch4_files = glob.glob('ch4/s*.fits')
    _, corrected_images = build_super_images(ch4_files)
    super_image_ch4 = corrected_images[0]

    # Now register all the super images to the J band image
    shifted_h = register_image(super_image_j, super_image_h)
    shifted_ch4 = register_image(super_image_j, super_image_ch4)

    # And write it out, with care to scoop some helpful header information
    hdu_h = fits.HDUList([fits.PrimaryHDU(shifted_h)])
    hdu_h[0].header = fits.open(h_files[0])[0].header

    hdu_j = fits.HDUList([fits.PrimaryHDU(super_image_j)])
    hdu_j[0].header = fits.open(j_files[0])[0].header

    hdu_ch4 = fits.HDUList([fits.PrimaryHDU(shifted_ch4)])
    hdu_ch4[0].header = fits.open(ch4_files[0])[0].header

    hdu_h.writeto('neptune_h.fits')
    hdu_j.writeto('neptune_j.fits')
    hdu_ch4.writeto('neptune_ch4.fits')
    

