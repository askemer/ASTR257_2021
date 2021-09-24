""" Not my best work I'm sorry I am very tired."""

## -- IMPORTS
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np

## -- FUNCTIONS

def build_wavelength_calibration(calibration_file):
    """ Creates a function to convert the wavelength calibration
    given a calibrationf ile."""

    wv_to_pix = build_wv_to_pix_dict(calibration_file)
    wv = list(wv_to_pix.keys())
    pix = np.array([wv_to_pix[key] for key in wv])

    b2, b1, m = np.polyfit(pix, wv, 2)

    def wv_to_pix(pix):
        return m + b1*pix + b2*pix**2
    
    return wv_to_pix
    

def build_wv_to_pix_dict(fits_file):
    """ Given we know the lines and where they fall we take the maxes 
    and return a dictionary of the pixel vs line values."""

    data = fits.getdata(fits_file)
    sums = np.sum(data, axis=0)

    wv_cal = {}
    wv_cal[5085.82] = np.where(sums == np.max(sums[1920:1940]))[0][0]
    wv_cal[4799.92] = np.where(sums == np.max(sums[1640:1680]))[0][0]
    wv_cal[4358.33] = np.where(sums == np.max(sums[1220:1260]))[0][0]
    wv_cal[4046.56] = np.where(sums == np.max(sums[920:940]))[0][0]
    wv_cal[3610.51] = np.where(sums == np.max(sums[450:500]))[0][0]
    wv_cal[3650.15] = np.where(sums == np.max(sums[500:550]))[0][0]
   
    return wv_cal


def plot_solar_lines(solar_wv, solar_sums):
    """ Plots up the solar spectrum and the lines all pretty."""

    lines = [3933.66, 3968.47, 4307.8, 4861.34]
    colors = ['cyan', 'teal', 'green', 'dodgerblue', 'skyblue']
    
    plt.clf()
    plt.plot(solar_wv, solar_sums, color='gray')
    
    for index, line in enumerate(lines):
        plt.axvline(line, color=colors[index], linestyle='--', label=str(line))

    plt.legend()
    plt.savefig('solar_spectrum_with_lines.png')


def sun_stuff(calibration_file, solar_file):
    """Main function to turn the crank."""

    wv_to_pix = build_wavelength_calibration(calibration_file)

    solar_data = fits.getdata(solar_file)
    solar_sums = np.sum(solar_data, axis=0)
    solar_wv = wv_to_pix(np.array([n for n in range(len(solar_sums))]))

    plot_solar_lines(solar_wv, solar_sums)

## -- RUN

if __name__ == "__main__":
    
    sun_stuff('./spectroscopy_data/arcs.fits', './spectroscopy_data/solar_spectrum.fits')



