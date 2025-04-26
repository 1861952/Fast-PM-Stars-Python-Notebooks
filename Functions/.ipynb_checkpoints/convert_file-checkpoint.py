import numpy as np
import pandas as pd

from astropy.io import fits
from astropy import units as u
from astropy.table import Table
from astropy import constants as const

#Scipy imports
from scipy.interpolate import InterpolatedUnivariateSpline, CubicSpline, Akima1DInterpolator
from scipy import integrate
from scipy.stats import powerlaw

from pathlib import Path
from Functions.f_nu_models import passband_flux, get_and_compare_mag


# these things don't change the value because they get cancelled for color (only depends on SED which depends on temperature)
R_sun = (1*u.Rsun).cgs.value
d = (10*u.pc).cgs.value
# F_0 in AB system, which is constant
F_0 = (3631.00*u.Jy).cgs.value

# Reads in Ascii file and returns wavelengths in cm and transmissions
def ascii_to_array(filename):
    # Read the text file into a DataFrame
    data = pd.read_csv(filename, sep=' ', header=None, names=['Wavelength', 'Transmission'])

    # wavelength in Angstrom
    column1 = (data['Wavelength'].values)*u.AA
    # Probability of photon going through
    column2 = data['Transmission'].values
    return column1, column2

# Displays .fits file as a table
def fits_to_table(fileURL):
    hdulist = fits.open(fileURL)
    data = Table(hdulist[1].data)
    print(data)

# Converts units of f_lambda in fits file to per cm, returns new converted table of valeus
def convert_fits_to_cm(fileURL):
    hdul = fits.open(fileURL)
    table_data = hdul[1].data.copy()
    log_g_list = table_data.dtype.names[1:]
    table = Table(table_data)
    table_data['WAVELENGTH'] *= 1e-8
    for log_g in log_g_list:
        table_data[log_g] *= 1e8   
    return table_data

# Reads in Fits file and returns actual values of wavelength and flux
def fits_to_wav_flux(fileURL, log_g):
    table_data = convert_fits_to_cm(fileURL)
    return table_data['WAVELENGTH']*u.cm, table_data[log_g]


# Reads in Fits file and returns spline for flux at surface (Phoenix Model)

def fits_to_spline(fileURL, log_g):
    hdul = fits.open(fileURL)
    data_fits = hdul[1].data
    # Wavelengths from fits table
    wavLen = data_fits['WAVELENGTH']
    # f_lamda values for specific log(g) 
    f_lambda = data_fits[log_g]
    # Makes spline based table
    # convert from Angstrom to cm
    wavelengthArr_cm = 1e-8*wavLen
    f_lambda_table_per_cm = f_lambda*1e8  
    # Returns spline at surface based on Phoenix model
    # convert from f_lambda to f_nu
    return InterpolatedUnivariateSpline(wavelengthArr_cm, f_lambda_table_per_cm*wavelengthArr_cm**2/const.c.cgs.value, k=1, ext="zeros")

def wav_trans_to_band_flux(fileURL, log_g, transmission_curve_file):
    wavelengths_band, transmission_band = ascii_to_array(transmission_curve_file)
    wavelengths_band_cm = wavelengths_band.cgs
    trans_band_spline = InterpolatedUnivariateSpline(wavelengths_band_cm, transmission_band, k=1, ext="zeros")
    wavs_actual, flux_actual = fits_to_wav_flux(fileURL, log_g)
    band_flux = passband_flux(flux_actual, trans_band_spline(wavs_actual.value), wavs_actual)
    return band_flux

def get_band_flux_dict(M_H, T_eff, log_g):
    band_flux_dict = {}
    bands = ["z", "Y", "i", "r"]
    for b in bands:
        trans_curve_band = f'CTIO/CTIO_DECam.{b}.txt'
        filename = f'All_Phoenix_Models/phoenixm{str(M_H).replace(".", "")}_{T_eff}.fits'
        flux = wav_trans_to_band_flux(filename, log_g, trans_curve_band)        
        band_flux_dict[b] = flux
    return band_flux_dict
    
# ----------------------------------------------------------
def fits_to_cubic_spline(fileURL, log_g):
    hdul = fits.open(fileURL)
    data_fits = hdul[1].data
    # Wavelengths from fits table
    wavLen = data_fits['WAVELENGTH']
    # f_lamda values for specific log(g) 
    f_lambda = data_fits[log_g]
    # Makes spline based table
    # convert from Angstrom to cm
    wavelengthArr_cm = 1e-8*wavLen
    f_lambda_table_per_cm = f_lambda*1e8  
    # Returns spline at surface based on Phoenix model
    # convert from f_lambda to f_nu
    return CubicSpline(wavelengthArr_cm, f_lambda_table_per_cm*wavelengthArr_cm**2/const.c.cgs.value)

def fits_to_akima_spline(fileURL, log_g):
    hdul = fits.open(fileURL)
    data_fits = hdul[1].data
    # Wavelengths from fits table
    wavLen = data_fits['WAVELENGTH']
    # f_lamda values for specific log(g) 
    f_lambda = data_fits[log_g]
    # Makes spline based table
    # convert from Angstrom to cm
    wavelengthArr_cm = 1e-8*wavLen
    f_lambda_table_per_cm = f_lambda*1e8  
    # Returns spline at surface based on Phoenix model
    # convert from f_lambda to f_nu
    return Akima1DInterpolator(wavelengthArr_cm, f_lambda_table_per_cm*wavelengthArr_cm**2/const.c.cgs.value)

# __________________________________________________________
# Finds nearest value (absolute value) in an array
def find_nearest_value(target, array):
    target = float(target)
    diff_arr = np.abs(np.asarray(array)-target)
    index = (np.abs(array - target)).argmin()
    return array[index]
# Creates a dictionary where key is (T_eff, log(g)) and value is the file name

def create_dict(foldername, glob_str):
    # Get folder
    path = Path(foldername)
    # path.glob(glob_str)

    # create dictionary : key is (T, log(g)) value is 'filename'
    file_dict = dict()
    for filepath in path.glob("sp_*"):
        filename = filepath.name
        fr = filepath.name.removeprefix("sp_t").removesuffix("nc_m-0.5.txt")
        T_eff = float(fr[0:4])
        g_mks = float(fr[5:])
        log_g_cgs = np.log10(g_mks*float(100))
        file_dict.update({(T_eff, log_g_cgs): filename})
        
    # np.unique automatically sorts
    T_eff_arr = np.unique([temp for temp, _lgg in file_dict])
    log_g_arr = np.unique([lgg for _temp, lgg in file_dict])
    T_eff_arr = np.asarray(T_eff_arr)
    log_g_arr = np.asarray(log_g_arr)
    return T_eff_arr, log_g_arr, file_dict


# Finds correct Sonora model file name based on T_eff and log(g) input
def find_file(T_eff_actual, log_g_actual, foldername, glob_str):
    T_eff_array, log_g_array, file_dict = create_dict(foldername, glob_str) 
    closest_T_eff = find_nearest_value(T_eff_actual, T_eff_array)
    closest_log_g = find_nearest_value(log_g_actual, log_g_array)
    g_mks = np.round(float(10)**(closest_log_g-float(2)))
    filename = f"{foldername}/sp_t{int(closest_T_eff)}g{int(g_mks)}nc_m-0.5.txt" 
    return filename

# Read in Sonora Model Spectra Files
def sonora_to_spline(filename):
    data = pd.read_csv(filename, skiprows=3, sep='\s+', names=['Wavelength', 'Flux'])
    # UNABLE TO DO ASTROPY UNITS FOR BOTH
    # convert from microns to centimeters
    wavelengths = data['Wavelength']*1e-4
    wavelengths = wavelengths[::-1]
    # observer flux
    Fluxes = data['Flux']
    Fluxes = Fluxes[::-1]
    # Makes spline based on table
    # Returns spline at surface based on Phoenix model
    # (already in f_nu so don't need to convert)
    return InterpolatedUnivariateSpline(wavelengths, Fluxes, k=1, ext="zeros")
