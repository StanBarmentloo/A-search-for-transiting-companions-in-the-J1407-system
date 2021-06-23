# -*- coding: utf-8 -*-
#Christiaan Dik and Stan Barmentloo 2020
#Bury Unwanted Trends to Correct and Horizontalize Every lightcuRve (BUTCHER)

import numpy as np
import scipy.signal as signal
from scipy.optimize import curve_fit
from astropy.timeseries import LombScargle
#import bootstrap
import matplotlib.pyplot as plt

def mag_to_flux(mag):
    """
    Converts a given array of magnitude values to normalized flux values.
    
    Parameters
    ----------
    mag : array-like
        The magnitude data.
    
    Returns
    -------
    flux : array-like
        The normalized flux data
    """
    flux = 10**(-0.4*mag)/np.mean(10**(-0.4*mag))
    return flux

def emag_to_eflux(mag, emag):
    """
    Converts a given array of uncertainty on magnitudes to normalized uncertainty on flux.
    
    Parameters
    ----------
    mag : array-like
        The magnitude data.
    emag : array-like
        The errors on the magnitude data. 
        Needs to have the same shape as mag.
        
    Returns
    -------
    eflux : array-like
        The error on the flux, normalized.
    """
    eflux = np.sqrt((-0.4*np.log(10)*np.exp(-0.4*np.log(10)*mag)*emag)**2)/np.mean(10**(-0.4*mag))
    return eflux


def long_correct(time, flux, eflux):
    """
    Fits a line to the data and removes this trend.
    
    Parameters
    ----------
    time : array-like
        The temporal data of the flux measurements.
    flux : array-like
        The flux data.
        Needs to have the same shape as time.
    eflux : array-like
        The error on the flux data.
        Need to have the same shape as time and flux.
        
    Returns
    -------
    corr_flux : array-like
        The flux data with a linear correction.
    """
    coef, cov = np.polyfit(time, flux, 1, w=1/(eflux**2), cov=True)
    corr_flux = flux/(coef[0]*time + coef[1])
    
    corr_flux/=np.mean(corr_flux)
    
    return corr_flux
    

def short_correct(time, flux, eflux, periods=np.linspace(3.1,3.3,500), gap_size=20, max_chunk_duration=500, return_durations=False, return_errors=False):
    """
    Divides the data into parts with gaps between them with no data (of a given gap size).
    The period with the strongest signal in the given period range is calculated.
    To the data folded over this period a sine is fitted and removed.
    Returns the total corrected flux and the fitted short periods.
    
    Parameters
    ----------
    time : array-like
        The temporal data of the flux measurements (in days)
    flux : array-like
        The flux data.
        Needs to have the same shape as time.
    eflux : array-like
        The error on the flux data.
        Need to have the same shape as time and flux.
    periods : array-like, optional
        The range of periods (in days) over which calculate the power in the periodogram.
        Default is 3 to 3.5 days with 500 intervals.
    gap_size: int or float, optional
        Size of the gap bewteen the data (in days) after which it is split.
        Default is 20 days.
    max_chunk_duration: int or float, optional
        Max duration of a piece (in days) before a new one is started.
        Default is 500 days.
    return_durations: Bool, optional
        Whether the timestamps of the chunks should be returned.
        Default is False.
    return_errors: Bool, optional
        Whether the Monte Carlo uncertainty of every period should be returned.
        Default is False.
        
    Returns
    -------
    final_flux :  array-like
        The flux corrected for the different short term periods and normalized.
    short_periods: array-like
        An array containing the best fitting short-term periods for the different pieces of data.
    durations: array-like, optional
        An array of begin and end times of the pieces of the data.
    errors: array-like, optional
        an array containing the uncertainties corresponding to short_periods.
        
    """
    angfreq = 2*np.pi/periods #Converts the periods to angular frequency
    short_corr_flux = np.copy(flux)

    short_periods = []
    FWHMs = []

    durations = []
    
    chunk_beginnings = [0]
    chunk_ends = []
    cutoffs = []

    for j in range(len(time)-1):
        #Look for gaps larger than gap_size days to identify 'chunks' of continuous measurements
        if (time[j]-time[j-1] > gap_size) or ((j-1) in cutoffs):
            chunk_beginnings.append(j)
        if time[j+1]-time[j] > gap_size:
            chunk_ends.append(j)
        #If the length of the chunk is greater than max_chunk_duration, start a new chunk
        if len(cutoffs) == 0 and len(chunk_ends) == 0:
            if time[j]-time[0] >= max_chunk_duration:
                chunk_ends.append(j)
                cutoffs.append(j)
        elif len(cutoffs) == 0 and len(chunk_ends) != 0:
            if time[j]-time[chunk_beginnings[-1]] >= max_chunk_duration:
                chunk_ends.append(j)
                cutoffs.append(j)
        else:
            if time[j]-time[cutoffs[-1]] >= max_chunk_duration:
                chunk_ends.append(j)
                cutoffs.append(j)
    chunk_ends.append(len(time)-1)

    found_amp = []  
    found_amp_err = []  
    found_period_err = [] 
    
    for i in range(len(chunk_beginnings)):
        if chunk_ends[i]-chunk_beginnings[i] > 6:
            chunk_time = time[chunk_beginnings[i]:chunk_ends[i]]
            chunk_flux = flux[chunk_beginnings[i]:chunk_ends[i]]
            chunk_eflux = eflux[chunk_beginnings[i]:chunk_ends[i]]

            frequency, power = LombScargle(chunk_time, chunk_flux-np.mean(chunk_flux)).autopower(minimum_frequency=1/(np.max(periods)), maximum_frequency=1/(np.min(periods)), samples_per_peak=10)
            
            
            
            found_period = (1/(frequency[np.argmax(power)]))
            periods = 1/np.array(frequency)
            short_period = float(found_period)

            if return_errors:
                found_err = monkey.main(chunk_time, chunk_flux, chunk_eflux, found_period)
                found_period_err.append(found_err)
            
            #Scipy lombscargle method below
            #pgram = signal.lombscargle(chunk_time, chunk_flux-np.mean(chunk_flux), angfreq, normalize=True)
            #short_period = periods[np.argmax(pgram)]

            #Determining the FWHM of the peak:
            #HM = 0.5*np.max(power)
            max_index = np.argmax(power)
            #left_HM = 0
            #right_HM = len(power)-1
            
            #Half max left:
            #for j in range(max_index):
            #    if (power[j+1] >= HM) and (power[j] < HM):
            #        left_HM = j+1

            #Half max right:
            #for j in range(len(power[max_index:(len(power)-1)])):
            #    if (power[max_index+j] >= HM) and (power[max_index+j+1] < HM):
            #        right_HM = max_index+j
            #        break

            #FWHM = 1/frequency[left_HM] - 1/frequency[right_HM]

            #plt.plot(frequency, power)
            #plt.axvline(frequency[left_HM])
            #plt.axvline(frequency[right_HM])
            #plt.axvline(frequency[max_index])
            #plt.show()
            
            short_periods.append(short_period)
            #FWHMs.append(FWHM)
            
            fold_time = chunk_time%short_period
            def sine(x, a, b, d):
            #Creates a sine curve with parameters a, b, d but with fixed c (period)
                return a+b*np.sin((2*np.pi/short_period)*(x-d))
        
            #Fit a sine curve with this set period
            popt, pcov = curve_fit(sine, fold_time, chunk_flux, sigma=chunk_eflux, absolute_sigma = True)
            found_amp.append(abs(popt[1]))
            found_amp_err.append(np.sqrt(abs(pcov[1][1])))
            chunk_flux = chunk_flux/sine(chunk_time, *popt)
            if len(chunk_flux) > 30: #minimal amount of points needed to fit sine
                short_corr_flux[chunk_beginnings[i]:chunk_ends[i]] /= sine(chunk_time, *popt) #remove the trend
    
    final_flux = short_corr_flux/np.mean(short_corr_flux)
    short_periods = np.array(short_periods)
    FWHMs = np.array(FWHMs)
    errors = np.array(found_period_err)
    
    if return_durations:
        for k in range(len(chunk_beginnings)):
            if chunk_ends[k]-chunk_beginnings[k] > 6:
                durations.append([time[chunk_beginnings[k]], time[chunk_ends[k]]])
        if return_errors:
            return final_flux, short_periods, durations, errors
        return final_flux, short_periods, durations

    if return_errors:
        return final_flux, short_periods, errors
    
    return final_flux, short_periods

def rand_short_period_finder(time, flux, eflux, sample_size=60, n_samples=100, gap_size=20, periods=np.linspace(3.1,3.3,500)):
    """
    Calculates the short 3.2 day rotation period at different times in the data.
    
    Parameters
    ----------
    time : array-like
        The temporal data of the flux measurements (in days)
    flux : array-like
        The flux data.
        Needs to have the same shape as time.
    eflux : array-like
        The error on the flux data.
        Need to have the same shape as time and flux.
    sample_size : int, optional
        The number of days over which a rotation period should be fitted.
        Default is 60 days.
    n_samples : int, optional
        The number of samples that should be taken from the data.
        Default is 100.
    periods : array-like, optional
        The range of periods (in days) over which calculate the power in the periodogram.
        Default is 3 to 3.5 days with 500 intervals.

    Returns
    -------
    midpoints : array-like
        An array containing the midpoints (in JD) over which the periods have been calculated.
    short_periods : array-like
        An array containing the calculated rotation periods.
    errors : array-like
        An array containing the errors on short_periods.
    """
    
    chunk_beginnings = [0]
    chunk_ends = []
    failed_points = []

    for j in range(len(time)-1):
        #Look for gaps larger than gap_size days to identify 'chunks' of continuous measurements
        if (time[j]-time[j-1] > gap_size) and (len(chunk_ends) != 0):
            chunk_beginnings.append(j)
        if time[j+1]-time[j] > gap_size:
            chunk_ends.append(j)
        
    chunk_ends.append(len(time)-1)

    
    midpoints, short_periods, found_period_err = [], [], []

    for k in range(len(chunk_beginnings)):
        if chunk_ends[k]-chunk_beginnings[k] > sample_size+1:
            left_edge = time[chunk_beginnings[k]]
            right_edge = time[:chunk_ends[k]][-1]-sample_size
            pickable_mask = (time < right_edge) * (time > left_edge)
            pick_time = time[pickable_mask]

            if len(pick_time) != 0:
                
                start_indices = np.arange(chunk_beginnings[k], chunk_ends[k]-sample_size, sample_size)
                
            else:
                start_indices = []
            for i in start_indices:
                end_i = np.argmin(abs(time-(time[i]+sample_size)))
                chunk_time = time[i:end_i]
                chunk_flux = flux[i:end_i]
                chunk_eflux = eflux[i:end_i]

                if (end_i-i)>10:
                    frequency, power = LombScargle(chunk_time, chunk_flux-np.mean(chunk_flux)).autopower(minimum_frequency=1/(np.max(periods)), maximum_frequency=1/(np.min(periods)), samples_per_peak=10)
                    found_period = (1/(frequency[np.argmax(power)]))
                    periods = 1/np.array(frequency)
                    short_period = float(found_period)
                    midpoint = np.median(time[i:end_i])
                    found_err = monkey.main(chunk_time, chunk_flux, chunk_eflux, found_period)
                    if (found_err != 'Fail'):
                        found_period_err.append(found_err)                
                        midpoints.append(midpoint)
                        short_periods.append(short_period)
    midpoints, short_periods, errors = np.array(midpoints), np.array(short_periods), np.array(found_period_err)
    return midpoints, short_periods, errors
