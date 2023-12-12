# Functions for calculating survey-related quantities

from __future__ import print_function

__author__ = "E.A.K. Adams"

"""
Collection of functions for calculating survey-related quantities
The focus is on being able to convert survey sensitivities to 
comparable quantities (e.g., uniform channel size)
"""

import numpy as np
import astropy.units as u
from astropy import constants as const

# globally define rest frequency
rest_freq = 1420.405752 * u.MHz

def get_nhi(sensitivity,linewidth,beam_maj,beam_min):
    """
    Get N_HI for given sensitivity
    Inputs:
    - sensitivity: Quantity, sensitivity matched to linewidth
    - linewidth: Quantity, linewidth for N_HI calc
    - beam_maj: Quantity, beam major axis
    - beam_min: Quantity, beam minor axis
    """
    # get beam area
    omega_beam = get_beam_area(beam_maj,beam_min)
    # convert to brightness temp
    tb = sensitivity.to(u.K,u.brightness_temperature(rest_freq,beam_area=omega_beam))
    # make sure linewidth is km/s (assume radio), not Hz units
    if 'Hz' in linewidth.unit.to_string():
        linewidth = convert_chan_freq_vel(linewidth)

    #convert to N_HI
    nhi = 1.823e18 * tb * linewidth  / (u.K * u.km / u.s) / (u.cm * u.cm)
    return nhi

def get_beam_area(beam_maj, beam_min):
    """
    Calculate beam area, provided FWHM of maj and min axis
    Inputs are quantities - have angular units
    """
    fwhm_to_sigma = 1. / (8 * np.log(2))**0.5
    maj_sigma = beam_maj * fwhm_to_sigma
    min_sigma = beam_min * fwhm_to_sigma
    beam_area = 2 * np.pi * min_sigma * maj_sigma
    
    return beam_area

def get_sens_for_res(sens,spec_res,desired_res):
    """
    Get the sensitivity for a desired spectral resolution
    Inputs:
    - sens: Quantity, sensitivey at given spec_res
    - spec_res: Quantity, spec. resolution for input sensitivity
    - desired_res: Quantity, desired spec res.
    Outputs:
    - new_sens: Quantity, sensitivity for desired_res
    """
    #first check that have same units for res
    #otherwise have to us an equivalency
    if not spec_res.unit == desired_res.unit:
        print("Units don't match")
        print("Checking if they are convertible")
        try:
            desired_res.to(spec_res.unit)
        except:
            print("Not convertible so change desired_res in freq/vel")
            desired_res = convert_chan_freq_vel(desired_res)
            # force to same units
            desired_res=desired_res.to(spec_res.unit)
            print(desired_res)

    #calculate the new sensitivity
    new_sens = sens * np.sqrt(spec_res / desired_res)
    print(new_sens)

    #return the new sensitivity value
    return new_sens
        
def convert_chan_freq_vel(chan_res):
    """
    Takes an input channel resolution, either in freq or velocity
    Will return matching resolution in other dimensions
    Presumes radio definition; can consider as a param for future
    Inputs:
    chan_res (Quantity): Input resolution, with units
    Outputs:
    new_res (Quantity): Output resolution, with units, opp of input
    """
    #define conversion
    freq_to_vel = u.doppler_radio(rest_freq)
    #check if in frequency; if so, get vel
    if 'Hz' in chan_res.unit.to_string():
        offsetfreq = rest_freq - chan_res
        vel = rest_freq.to(u.km/u.s, equivalencies = freq_to_vel)
        offsetvel = offsetfreq.to(u.km/u.s, equivalencies = freq_to_vel)
        new_res = offsetvel-vel

    else:
        print("Not frequency unit, presuming velocity")
        offsetfreq = chan_res.to(u.kHz, equivalencies = freq_to_vel)
        new_res = rest_freq-offsetfreq

    return new_res

def convert_vel_freq(quantity):
    """
    Take an input quantity that is either vel or freq and return the other
    Default for HI line
    Use radio convention (could set as param)
    """
    freq_to_vel = u.doppler_radio(rest_freq)
    if 'Hz' in quantity.unit.to_string():
        new_quantity = quantity.to(u.km/u.s, equivalencies = freq_to_vel)
    else:
        #presume in vel, use try/except in case
        try:
            new_quantity = quantity.to(u.MHz,  equivalencies = freq_to_vel)
        except:
            print("Units of input not recognized")
            new_quantity = np.nan

    return new_quantity
