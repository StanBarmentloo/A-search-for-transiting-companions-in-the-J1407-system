#!/usr/bin/env python
# coding: utf-8

# In[8]:


import numpy as np
import bootstrap
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.optimize import curve_fit
from astropy.timeseries import LombScargle


# In[ ]:

def sine(x, period, amplitude, offset, base):
    return amplitude*np.sin((2*np.pi/period)*(x-offset))+base



# In[ ]:


def alt_period_finder(time, flux, eflux, gap_size = 20, sample_size = 30, absolute_minimal_chunk_length =20, periods = np.linspace(3.1,3.3,500), points_per_cycle = 2, plot = False):
    
    chunk_beginnings = [time[0]]
    chunk_ends = []
    failed_points = []

    for j in range(len(time)-1):
        #Look for gaps larger than gap_size days to identify 'chunks' of continuous measurements
        if (time[j]-time[j-1] > gap_size) and (len(chunk_ends) != 0):
            chunk_beginnings.append(time[j])
        if time[j+1]-time[j] > gap_size:
            chunk_ends.append(time[j])
        
    chunk_ends.append(time[-1])

    
    midpoints, short_periods, found_period_err = [], [], []
    
    
    max_powers = []
    
    for k in range(len(chunk_beginnings)):
        chunk_length = chunk_ends[k]-chunk_beginnings[k]
        if chunk_length > absolute_minimal_chunk_length:
            for l in range(np.max((int(chunk_length/sample_size), 1))):
                #print("Hier is l ", l)
                segment_mask = (chunk_beginnings[k]+l*sample_size < time) * (time < chunk_beginnings[k]+(l+1)*sample_size)
                
                segment_time = time[segment_mask]
                segment_flux = flux[segment_mask]
                segment_eflux = eflux[segment_mask]
                    
                if len(segment_time) > (float(sample_size)/3.2)*points_per_cycle: 
                    frequencies = 1/np.linspace(2, 5, 3000)
                    power = LombScargle(segment_time, segment_flux-np.mean(segment_flux), dy = segment_eflux).power(frequencies)
                
                    if 3.1 < 1/frequencies[np.argmax(power)] < 3.3:
                    
                    
                    
                        frequency, power = LombScargle(segment_time, segment_flux-np.mean(segment_flux)).autopower(minimum_frequency=1/(np.max(periods)), maximum_frequency=1/(np.min(periods)), samples_per_peak=10)
                        #popt, pcov = curve_fit(sine, segment_time, segment_flux, p0=[3.2, 0.05, 13, 1.0], bounds=((3.1, 0, 0, 0.8), (3.3, 0.2, np.inf, 1.2)), absolute_sigma = True, sigma=segment_eflux, max_nfev = 2000) #sigma = 1/(segment_eflux**2)
                        #print(popt)
                        #time_space = np.linspace(np.min(segment_time), np.max(segment_time), 1000)
                        #plt.plot(time_space, sine(time_space, *popt))
                        #plt.scatter(segment_time, segment_flux)
                        #plt.xlim(np.min(segment_time), np.max(segment_time))
                        #plt.ylim(0.8, 1.2)
                        #plt.show()

                        found_period = float((1/(frequency[np.argmax(power)])))
                        midpoint = np.median(segment_time)

                        max_powers.append(power[np.argmax(power)])

                        if found_period < 3.17:
                            plt.show()
                            print(found_period)
                            print('de max power is', '{0:.2f}'.format(power[np.argmax(power)]))
                            plt.plot(1/frequency, power)
                            plt.show()

                            plt.scatter(segment_time, segment_flux)
                            plt.plot(segment_time, segment_flux, c = 'r')
                            plt.show()

                        found_err = bootstrap.main(segment_time, segment_flux, segment_eflux, found_period)
                        if (found_err != 'Fail'):
                            found_period_err.append(found_err)                
                            midpoints.append(midpoint)
                            short_periods.append(found_period)

                        #-------------------------------------------------------------------------------
                        if plot == True and k == 7 and l == 0:
                            plt.show()
                            fig, ax = plt.subplots(2, 1,  gridspec_kw={'height_ratios': [3, 1.5]})

                            gs1 = gridspec.GridSpec(2, 1)
                            gs1.update(wspace=0.025, hspace=0.05)

                            ax[0].errorbar(segment_time, segment_flux, segment_eflux, fmt = '.')
                            ax[0].set_xlabel('MJD (Days)')
                            ax[0].set_ylabel('Flux (Normalised Flux)')



                            ax[1].plot(1/frequency, power)
                            ax[1].axvline(x = found_period, c= 'r', linestyle = '--', alpha = 0.5)
                            ax[1].set_xlabel('Periodicity (Days)')
                            ax[1].set_ylabel('Signal Power (Arbitrary Units)')
                            plt.tight_layout()
                            fig = plt.gcf()
                            fig.set_size_inches(12, 9)
                            plt.show()


                        #-------------------------------------------------------------------------------
                        if 3.11 < found_period < 3.29: 
                            try:
                                period_format = '{0:.' +str(len(str('{0:.1g}'.format(found_err)))-1)+ 'g}'
                            except:
                                period_format = '{0:.1g}'
                            try:
                                print(str('{0:.3f}'.format(found_period)) +' & '+ str('{0:.3f}'.format(found_err)) +' & '+ str(int(np.min(segment_time)))+ ' & '+str(int(np.max(segment_time))) + ' \\\\')
                            except: 
                                pass
                        #-----------------------------------------------------------------------------------

                elif len(segment_time)/(segment_time[-1]-segment_time[0]) > points_per_cycle :
                    frequencies = 1/np.linspace(2, 5, 3000)
                    power = LombScargle(segment_time, segment_flux-np.mean(segment_flux), dy = segment_eflux).power(frequencies)
                
                    if 3.1 < 1/frequencies[np.argmax(power)] < 3.3:
                    
                    
                        frequency, power = LombScargle(segment_time, segment_flux-np.mean(segment_flux)).autopower(minimum_frequency=1/(np.max(periods)), maximum_frequency=1/(np.min(periods)), samples_per_peak=10)
                        #try:
                        #    popt, pcov = curve_fit(sine, segment_time, segment_flux, p0=[3.2, 0.05, 13, 1.0], bounds=((3.1, 0, 0, 0.8), (3.3, 0.2, np.inf, 1.2)), absolute_sigma = True, sigma=segment_eflux, max_nfev = 2000) 
                        #    print(popt)
                        #    time_space = np.linspace(np.min(segment_time), np.max(segment_time), 1000)
                        #    plt.plot(time_space, sine(time_space, *popt))
                        #    plt.scatter(segment_time, segment_flux)
                        #    plt.xlim(np.min(segment_time), np.max(segment_time))
                        #    plt.ylim(0.8, 1.2)
                        #    plt.show()
                        #except:
                        #    pass

                        found_period = float((1/(frequency[np.argmax(power)])))
                        midpoint = np.median(segment_time)

                        if found_period < 3.17:
                            plt.show()
                            print(found_period, 'het zit in de else')
                            print('de max power is', '{0:.2f}'.format(power[np.argmax(power)]))
                            plt.plot(1/frequency, power)
                            plt.show()

                            plt.scatter(segment_time, segment_flux)
                            plt.plot(segment_time, segment_flux, c = 'r')
                            plt.show()




                        found_err = bootstrap.main(segment_time, segment_flux, segment_eflux, found_period, tries = 300)
                        if (found_err != 'Fail'):
                            found_period_err.append(found_err)                
                            midpoints.append(midpoint)
                            short_periods.append(found_period)

                    #print("This segment was from day ", int(np.min(segment_time)), "to ", int(np.max(segment_time)))
                    #print("The found period was ", found_period)
                    #print("The found error was ", found_err)
                        if 3.11 < found_period < 3.29: 
                            try:
                                period_format = '{0:.' +str(len(str('{0:.1g}'.format(found_err)))-1)+ 'g}'
                            except:
                                period_format = '{0:.1g}'
                            try:
                                print(str('{0:.3f}'.format(found_period)) +' & '+ str('{0:.3f}'.format(found_err)) +' & '+ str(int(np.min(segment_time)))+ ' & '+str(int(np.max(segment_time))) + ' \\\\')
                            except: 
                                pass
                else:
                    print('{0:.2f}'.format(segment_time[0]), ' days.' )
                    #plt.plot(1/frequencies, power)
                    #plt.show()
                    
    #plt.show()        
    #plt.hist(max_powers, density = True, color= 'green')
    #plt.show()
    midpoints, short_periods, found_period_err = np.array(midpoints), np.array(short_periods), np.array(found_period_err)
    return midpoints, short_periods, found_period_err

