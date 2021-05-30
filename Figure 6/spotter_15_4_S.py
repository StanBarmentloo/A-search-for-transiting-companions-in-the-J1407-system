#!/usr/bin/env python
# coding: utf-8
#Searching PeriodOgrams for Transits To EmeRge (SPOTTER)


# In[1]:


import numpy as np
import matplotlib.pyplot as plt
#import batman
#import bls
import scipy.signal as signal
from astropy import units as u
from astropy.io import ascii
from astropy.timeseries import LombScargle
#from astropy.timeseries import BoxLeastSquares
#from transitleastsquares import transitleastsquares

r_jupiter = u.astrophys.jupiterRad.to(u.m)
r_star = 0.96*u.astrophys.R_sun.to(u.m)
m_star = 0.95*u.astrophys.M_sun.to(u.kg)
day_in_sec = u.day.to(u.second)
G = 6.6743*10**-11
resonances = [1./5.,1./4.,1./3., 1./2., 1., 2., 3., 4., 5.]


# In[2]:
def transit_time(period):
    a = float(((((period*day_in_sec)**2)*G*m_star)/(4*np.pi**2))**(1./3.))/r_star 
    return period/(np.pi*a)




#Generating noise and fake data
def noise_generator(times, noise_type, std):
    if noise_type == 'Gauss':
        return np.random.randn(times.size)*std #Generating Gaussian noise with a sigma of choice
    elif noise_type == 'Poisson':
        pass
        
def fake_data_generator(times, noise):
    data = np.ones_like(times) 
    data += noise
    return data


# In[3]:


#Generating a fake transit
def transit_generator(times, r_planet, t0, per, m_star, r_star, inc, ecc, w):
    params = batman.TransitParams()              #object to store transit parameters
    params.t0 = float(t0)                        #time of inferior conjunction
    params.per = float(per)                      #orbital period
    params.rp = float(r_planet*r_jupiter)/r_star #planet radius (in units of stellar radii)
    params.a = float(((((per*day_in_sec)**2)*G*m_star)/(4*np.pi**2))**(1./3.))/r_star   #semi-major axis (in units of stellar radii)
    params.inc = float(inc)                      #orbital inclination (in degrees)
    params.ecc = float(ecc)                      #eccentricity
    params.w = float(w)                          #longitude of periastron (in degrees)
    params.limb_dark = "nonlinear"               #limb darkening model "nonlinear"
    params.u = [1.2, -0.47, -0.22, 0.24]         #limb darkening coefficients [u1, u2, u3, u4] for T4500 log4.0, standard = [0.5, 0.1, 0.1, -0.1] 
    #print("The transit takes", params.per/(np.pi*params.a), "Days")


    t = times #times at which to calculate light curve 
    m = batman.TransitModel(params, t)    #initializes model
    #plt.scatter(t%per, m.light_curve(params), s=2)
    #plt.show()
    return m.light_curve(params)


# In[4]:


#Finding the transit using boxleastsquares algorithm
def transit_finder(times, transit_flux, yerr, min_p, max_p, plot, stats): #astropy
    model = BoxLeastSquares(np.array(times), transit_flux, dy = yerr)
    #periodogram = model.autopower(1, minimum_period= float(min_p), maximum_period= float(max_p))
    print (transit_time(min_p), transit_time(max_p))
    periodogram = model.power(np.linspace(min_p, max_p, 1000), np.linspace(transit_time(min_p), transit_time(max_p), 100))
    
    if plot == False: #If False, skip all plotting to save computing time
        pass
    else:
        plt.errorbar(times, transit_flux, yerr=yerr, fmt='.', markersize= 3, elinewidth=0.5)  #Plot the fake transit data
        plt.xlabel("Julian Date (Days)") 
        plt.ylabel("Normalised Flux")
        figure = plt.gcf()
        figure.set_size_inches(18, 10)
        plt.show()
        
        plt.plot(periodogram.period, periodogram.power) #Plot the periodogram with vertical lines at resonance peaks
        for i in resonances:
            if i*periodogram.period[np.argmax(periodogram.power)] > min_p and i*periodogram.period[np.argmax(periodogram.power)] < max_p:
                plt.axvline(i*periodogram.period[np.argmax(periodogram.power)], c='r', linestyle='--', alpha =0.5)
        plt.xlabel("Period (days)")
        plt.ylabel("Power")
        figure = plt.gcf()
        figure.set_size_inches(16, 8)
        plt.show()
        
    if stats == False: #If False, skip printing statistics to save computing time
        pass
    else:
        max_power_index = np.argmax(periodogram.power)
        print(model.compute_stats(periodogram.period[max_power_index], periodogram.duration[max_power_index], periodogram.transit_time[max_power_index]))  #erroring somewhere

    #print("Orbital period of", str(periodogram.period[np.argmax(periodogram.power)])[:5], "days")
    return periodogram.period[np.argmax(periodogram.power)]

# In[6]:
def alt_transit_finder(times, transit_flux, yerr, min_p, max_p, plot, stats): #kovacs
    nfreqs=int((1/float(min_p)-1/float(max_p))/0.0001)
    fmin = 1/float(max_p)
    results = bls.eebls(np.array(times), 
                        np.array(transit_flux),
                        u=np.zeros_like(times), 
                        v=np.zeros_like(times), 
                        nf=nfreqs,
                        fmin=1/float(max_p),
                        df=0.0001,
                        nb=1750, 
                        qmi=0.005, 
                        qma=0.12)
    
    if plot == False: #If False, skip all plotting to save computing time
        pass
    else:
        plt.errorbar(times, transit_flux, yerr=yerr, fmt='.', markersize= 3, elinewidth=0.5) #Plot the fake transit data
        plt.xlabel("Julian Date (Days)")
        plt.ylabel("Normalised Flux")
        figure = plt.gcf()
        figure.set_size_inches(18, 10)
        plt.show()
        
        freq_range = np.linspace(fmin, fmin+nfreqs*0.0001, nfreqs)
        plt.plot(1/freq_range, results[0]) #Plot periodogram with vertical lines at resonance peaks                           
        for i in resonances:
            if i*(1/freq_range)[np.argmax(results[0])] > min_p and i*(1/freq_range)[np.argmax(results[0])] < max_p:
                plt.axvline(i*(1./freq_range)[np.argmax(results[0])], c='r', linestyle='--', alpha =0.5)
        plt.xlabel("Period (days)")
        plt.ylabel("Power")
        figure = plt.gcf()
        figure.set_size_inches(16, 8)
        plt.show()

    #print("Orbital period of", str(results[1])[:5], "days")
    return results[1]

def tls(times, transit_flux, yerr, min_p, max_p, plot, stats): #tls 
    model = transitleastsquares(times, transit_flux, yerr)
    results = model.power(R_star = 0.96, R_star_min = 0.81, R_star_max = 1.11, M_star = 0.95, M_star_min = 0.85, M_star_max = 1.05, period_min = min_p, period_max = max_p, oversampling_factor = 1, use_threads = 36, show_progress_bar = True)
    #R_star = 0.99, R_star_min = 0.90, R_star_max = 1.10, M_star = 0.90, M_star_min = 0.8, M_star_max = 1.0,
    
    if stats == False:
        pass
    else:
        print('Period', format(results.period, '.5f'), 'd')
        #print(len(results.transit_times), 'transit times in time series:', \
        #['{0:0.5f}'.format(i) for i in results.transit_times])
        print('Transit depth', format(results.depth, '.5f'))
        print('Best duration (days)', format(results.duration, '.5f'))
        print('Signal detection efficiency (SDE):', results.SDE)
        
    if plot == False:
        pass
    else:
        plt.figure()
        ax = plt.gca()
        ax.axvline(results.period, alpha=0.4, lw=3)
        for n in range(2, 10):
            ax.axvline(n*results.period, alpha=0.4, lw=1, linestyle="dashed")
            ax.axvline(results.period / n, alpha=0.4, lw=1, linestyle="dashed")
        plt.ylabel(r'SDE')
        plt.xlabel('Period (days)')
        plt.title('Periodogram Final Data')
        plt.plot(results.periods, results.power, color='black', lw=0.5)
        plt.xlim(np.min(results.periods), np.max(results.periods))
        #plt.savefig('Periodogram_7min_talk.png', dpi = 300)

        plt.figure()
        plt.plot(results.model_folded_phase, results.model_folded_model, color='red')
        plt.errorbar(results.folded_phase, results.folded_y, fmt='.', color='blue', alpha=0.5, zorder=2) #Might error here
        plt.xlim(0.48, 0.52)
        plt.xlabel('Phase')
        plt.ylabel('Relative flux');
    
    #Picking out the three most likely periods, excluding possible insignificantly different values
    picker = 1
    top_4_SDE = [0]
    while len(top_4_SDE) < 4:
        if len(top_4_SDE) == 3:
            if abs((results.periods[list(results.power).index(sorted(list(results.power))[-picker])]-results.periods[list(results.power).index(top_4_SDE[-1])])/(results.periods[list(results.power).index(top_4_SDE[-1])])) < 0.05 or abs((results.periods[list(results.power).index(sorted(list(results.power))[-picker])]-results.periods[list(results.power).index(top_4_SDE[-2])])/(results.periods[list(results.power).index(top_4_SDE[-2])])) < 0.05:
                pass
            else:
                top_4_SDE.append(sorted(list(results.power))[-picker])
        elif len(top_4_SDE) == 1:
            top_4_SDE.append(sorted(list(results.power))[-1])
        else:
            if abs((results.periods[list(results.power).index(sorted(list(results.power))[-picker])]-results.periods[list(results.power).index(top_4_SDE[-1])])/(results.periods[list(results.power).index(top_4_SDE[-1])])) < 0.05:
                pass
            else:
                top_4_SDE.append(sorted(list(results.power))[-picker])
        picker += 1
        
    top_3_SDE = top_4_SDE[1:]
    top_3_periods = [results.periods[list(results.power).index(top_3_SDE[0])], results.periods[list(results.power).index(top_3_SDE[1])], results.periods[list(results.power).index(top_3_SDE[2])]]
    
    print(top_3_periods, top_3_SDE)
    
    return top_3_periods, top_3_SDE, results



#Main function running the entire series
def main(times, transit_flux, yerr, min_p, max_p, plot, stats, add_transit, transit_params, version):
    if add_transit == False:
        pass
    else:
        transit_flux += (transit_generator(times, *transit_params)-1) #Building a dataset including a fake transit
    if version == 'kovacs':
        results = alt_transit_finder(times, transit_flux, yerr, min_p, max_p, plot, stats)         
    elif version == 'astropy':
        results = transit_finder(times, transit_flux, yerr, min_p, max_p, plot, stats) #Searching for a transit in the inserted data set
    elif version == 'tls': 
        results = tls(times, transit_flux, yerr, min_p, max_p, plot, stats)
    return results

# In[ ]:




