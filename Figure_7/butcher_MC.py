# -*- coding: utf-8 -*-
#Christiaan Dik and Stan Barmentloo 2020
#Bury Unwanted Trends to Correct and Horizontalize Every lightcuRve (BUTCHER)

import numpy as np
import scipy.signal as signal
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from astropy.timeseries import LombScargle

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
    #plt.errorbar(time, flux, yerr=eflux, fmt='.', markersize= 3, elinewidth=0.5, label='Uncorrected Flux')
    #plt.plot(time, coef[0]*time + coef[1], linestyle = '--', c='red', label='Linear Fit')
    #plt.legend()
    #figure = plt.gcf()
    #figure.set_size_inches(18, 10)
    #plt.xlabel('MJD (days)')
    #plt.ylabel('Normalised Flux')
    #plt.title('Linear Fit to AAVSO Flux Data')
    #plt.savefig('Draft_Figure_1.png')
    #plt.show()
    
    corr_flux = flux/(coef[0]*time + coef[1])
    
    return corr_flux
    

def short_correct(time, flux, eflux, periods=np.linspace(3,3.5,500), gap_size=20, max_chunk_duration=500, return_durations=False):
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
        Default is 2 to 10 days with 1000 intervals.
    gap_size: int or float, optional
        Size of the gap bewteen the data (in days) after which it is split.
        Default is 20 days.
    max_chunk_duration: int or float, optional
        Max duration of a piece (in days) before a new one is started.
        Default is 500 days.
    return_durations: Bool, optional
        Whether the timestamps of the chunks should be returned.
        Default is False.
        
    Returns
    -------
    final_flux :  array-like
        The flux corrected for the different short term periods and normalized.
    short_periods: array-like
        An array containing the best fitting short-term periods for the different pieces of data.
    durations: array-like, optional
        An array of begin and end times of the pieces of the data.
        
    """
    angfreq = 2*np.pi/periods #Converts the periods to angular frequency
    short_corr_flux = np.copy(flux)

    short_periods = []

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
    
    found_amp = []      #This is different
    found_amp_err = []  #This is different
    found_period_err = [] #This is different 1
    
    for i in range(len(chunk_beginnings)):
        if chunk_ends[i]-chunk_beginnings[i] > 6:
            chunk_time = time[chunk_beginnings[i]:chunk_ends[i]]
            chunk_flux = flux[chunk_beginnings[i]:chunk_ends[i]]
            chunk_eflux = eflux[chunk_beginnings[i]:chunk_ends[i]]

            frequency, power = LombScargle(chunk_time, chunk_flux-np.mean(chunk_flux)).autopower(minimum_frequency=1/(np.max(periods)), maximum_frequency=1/(np.min(periods)))
            found_period = (1/(frequency[np.argmax(power)]))
            short_period = float(found_period)
            
            found_err = 0
            found_period_err.append(found_err) #This is different 1
            
            #Scipy lombscargle method below
            #pgram = signal.lombscargle(chunk_time, chunk_flux-np.mean(chunk_flux), angfreq, normalize=True)
            #short_period = periods[np.argmax(pgram)]

            short_periods.append(short_period)
            
            fold_time = chunk_time%short_period
            def sine(x, a, b, d):
            #Creates a sine curve with parameters a, b, d but with fixed c (period)
                return a+b*np.sin((2*np.pi/short_period)*(x-d))
        
            #Fit a sine curve with this set period
            popt, pcov = curve_fit(sine, fold_time, chunk_flux, p0 = [1, 0.02, 0], sigma=chunk_eflux, absolute_sigma = True)
            found_amp.append(abs(popt[1])) #This is different
            found_amp_err.append(np.sqrt(abs(pcov[1][1]))) #This is different
            chunk_flux = chunk_flux/sine(chunk_time, *popt)
            if len(chunk_flux) > 30: #minimal amount of points needed to fit sine
                short_corr_flux[chunk_beginnings[i]:chunk_ends[i]] /= sine(chunk_time, *popt) #remove the trend
    final_flux = short_corr_flux/np.mean(short_corr_flux)
    short_periods = np.array(short_periods)
    found_period_err = np.array(found_period_err)
    
    print("Average found amp =", sum(found_amp)/len(found_amp), len(found_amp)) #This is different
    print("Average found amp_err =", sum(found_amp_err)/len(found_amp_err), len(found_amp_err)) #This is different
    print("Found period errors are", found_period_err) #This is different 1
    if return_durations:
        for k in range(len(chunk_beginnings)):
            if chunk_ends[k]-chunk_beginnings[k] > 6:
                durations.append([time[chunk_beginnings[k]], time[chunk_ends[k]]])
        return final_flux, short_periods, durations, found_period_err #This is different 1
    
    return final_flux, short_periods

def long_sine_correct(time, flux, eflux, plot):
    def sine(x, a, b, c, d):
            #Creates a sine curve 
                return a+b*np.sin((2*np.pi/c)*(x-d))
    popt, pcov = curve_fit(sine, time, flux, sigma=eflux, absolute_sigma = True)
    print(popt[2])
    if plot == True:
        plt.errorbar(time, flux, eflux, fmt='.', markersize= 3, elinewidth=0.5)
        plt.plot(time, sine(time, *popt))
        figure = plt.gcf()
        figure.set_size_inches(18, 10)
        plt.show()
    return flux/(sine(time, *popt))
