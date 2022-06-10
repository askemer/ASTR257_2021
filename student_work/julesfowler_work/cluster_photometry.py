## -- IMPORTS
from astropy.io imports fits
from astropy.stats import mad_std
import matplotlib.pyplot as plt
import numpy as np
from photutils.aperture import aperture_photometry, CircularAperture
from photuils import centroid_1dg, centroid_2dg


## -- FUNCTIONS

# - do photometry on landolt fields
# - calculate zero points 
# - find sources in clusters
# - apply offset between scenes
# - do photometry in clusters
# - calculate magnitudes 
# - plot HR diagram

def build_hr_diagram():
    
    # build dictionary 
    landolt_dict = build_landolt_dictionary()
    
    # do cluster photometry 
    v_phot, b_phot = perform_cluster_photometry()

    # calculate magnitudes 
    v_mags = calculate_magnitude(v_phot, landolt_dict['V']['average_zp'])
    b_mags = calculate_magnitude(b_phot, landolt_dict['B']['average_zp'])

    # plot HR diagram
    plot_HR_diagram(v_mags, b_mags)

def perform_cluster_photometry(file_path):
    """ This is mostly a sketch, to be honest."""

    dat_B = fits.getdata(file_path + 'final_ngc6819_B.fits')
    dat_V = fits.getdata(file_path + 'final_ngc6819_V.fits')

    bkg_sigma_B = mad_std(dat_B)
    bkg_sigma_V = mad_std(dat_V)

    daofind_b = DAOStarFinder(fwhm=4., threshold=3. * bkg_sigma_B)
    daofind_v = DAOStarFinder(fwhm=4., threshold=3. * bkg_sigma_V)

    sources_b = daofind_b(dat_B)
    sources_v = daofind_v(dat_V)
    
    # this is rough I'm so sorry
    match_xv, match_yv, match_xb, match_yb = [], [], [], []
    for index, source in enumerate(sources_v):
        diff = sources_B['xcentroid'] - source['xcentroid']
        locs = np.where(np.abs(diff) < 1)[0]
        for val in locs:
            # apply cut for bad pixel column
            if np.abs(sources_B['ycentroid'][val] - source['ycentroid']) and source['xcentroid'] < 1025:
                match_xv.append(source['xcentroid'])
                match_yv.append(source['ycentroid'])
                match_xb.append(sources_B['xcentroid'][val])
                match_yb.append(sources_B['ycentroid'][val])
    
    # this gives you some doubles so take a set of the tuples
    tups_v = list(set([(match_xv[n], match_yv[n]) for n in range(len(match_xv))]))
    tups_b = list(set([(match_xb[n], match_yb[n]) for n in range(len(match_xb))]))

    apertures_v = CircularAperture(tups_v, r=4.)
    apertures_b = CircularAperture(tups_b, r=4.)

    phot_v = aperture_photometry(dat_B, apertures_v)
    phot_b = aperture_photometry(dat_V, apertures_b)

    # and another round of crossmatching 
    phot_vals_v, phot_vals_b, x, y = [], [], [], []
    for index, source in enumerate(phot_v):
        print(source)
        diff = np.array(phot_b['xcenter'] - source['xcenter'])
        locs = np.where(np.abs(diff) < 1)[0]
        for val in locs:
            if np.abs(phot_b['ycenter'][val] - source['ycenter']):
                phot_vals_v.append(source['aperture_sum'])
                phot_vals_b.append(phot_b['aperture_sum'][val])
    
    return np.array(phot_vals_v), np.array(phot_vals_b)



def build_landolt_dictionary(file_path)

    # calculate landolt radius
    landolt_radius = (7*(6.3/1024))/2
    
    landolt_B = fits.getdata(file_path + 'final_landolt_B.fits')
    landolt_V = fits.getdata(file_path + 'final_landolt_V.fits')

    # given known regions where landolt stars hang and their magnitudes
    landolt_dict = {}
    landolts = [('1965', {'V':11.419, 'B':11.419+1.710}, 
                         {'V':V_data[360:400, 320:370],
                          'B':B_data[360:400, 320:370]}), 
                ('1969', {'V':10.382, 'B':10.382+1.959}, 
                         {'V':V_data[520:560, 250:290],
                          'B':B_data[520:560, 250:290],}), 
                ('1925', {'V':12.388, 'B':12.388+0.395},
                         {'V':V_data[660:720, 830:880], 
                          'B':B_data[660:720, 830:880]})]
    bands = ['V', 'B']
    regions = ???
    for band in bands:
        avg = 0
        for landolt in landolts:
        
            target, magnitude, region = landolt
            landolt_dict[band] = {}
            
            landolt_dict[band][target] = {}
            fits_file = f'final_landolt_{band}.fits'
            data = fits.getdata(fits_file)
            
            landolt_dict[band][target]['region'] = region
            photometry = annulus_photometry(region, landolt_radius)
            if photometry < 1000:
                photometry = annulus_photometry(region, landolt_radius, centroid='centroid_2dg')
            landolt_dict[band][target]['photometry'] = photometry
            zp = calculate_zero_point(photometry, magnitude)
            landolt_dict[band][target]['zp'] = zp
            avg += zp
        landolt_dict[band]['average_zp'] = avg/3

    return landolt_dict


def annulus_photometry(region, radius, centroid='centroid_1dg'):
    
    centroid_method = centroid_1dg if centroid=='centroid_1dg' else centroid_2dg
    x, y = centroid_method(region)

    source_aperture = CircularAperture([(x,y)], r=radius)
    ring_aperture = CircularAnnulus([(x,y)], r_in=radius+radius/2, r_out=2*radius+radius/2)
    
    region_ratio = radius**2/((2*radius+radius/2)**2 - (radius+radius/2)**2)

    source_photometry = aperture_photometry(region, source_aperture)['aperture_sum'][0]
    ring_photometry = aperture_photometry(region, ring_aperture)['aperture_sum'][0]*region_ratio
    
    background_subtracted_photometry = source_photometry - ring_photometry
    
    return background_subtracted_photometry 


def calculate_magnitude(counts, zero_point):

    magnitude = -2.5*np.log10(counts) + zero_point
    return magnitude

def calculate_zero_point(counts, magnitude):

    zero_point = magnitude + 2.5*np.log10(counts)
    return zero_point

def plot_HR_diagram(v_mags, b_mags, mode='B'):

    plt.clf()
    mags = b_mags if mode=='B' else v_mags
    plt.scatter(b_mags-v_mags, mags, s=10, edgecolor='gray', facecolor='cyan',
            alpha=.5, s=20)
    plt.xlabel('B-V color [mag]')
    plt.ylabel(f'{mode} magnitude [mag]')
    plt.ylim(np.max(mags), np.min(mags))
    plt.savefig('hr_diagram.png')
