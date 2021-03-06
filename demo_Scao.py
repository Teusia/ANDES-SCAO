# -*- coding: utf-8 -*-
"""
Created on Fri Jan 28 11:45:35 2022

@author: anne.cheffot
@author: EstherSoria
"""
import numpy as np
import poppy as po
import statistics
import datetime
from pupil_ANDES import *
from sklearn.metrics import mean_squared_error
import math
import time
from astropy import units as u

def SCAOSim(wavelength,nPix,nScreen,SDIslandEffect,SDSegmentJitter,SD_WS
            ,mag,coro,pupilFileName='Tel-Pupil.fits', atmFilePrefix='screen'):
    '''
    Parameters
    ----------
    wavelength : float
        wavelength of the final image in meters
    nPix : int
        number of pixels across one image
    nScreen : int
        number of time light is propagated through an aberated screen.
    SDIslandEffect : float
        SD error of the island effect in meters
    SDSegmentJitter : float
        SD of the segment jitter in meters
    rmsWindshake : float
        root mean square of the tip tilt induced by the winshake in meters
    seed : int
        seed for the randome number generator
    atmWorseningFactor : float
        a number to change the atmospheric residuals amplitude (no unit)
    pupilFileName : string, optional
        The name of the file in which the pupil shape is stored. The default is 
        'Tel-Pupil.fits'.
    atmFilePrefix : string, optional
        The atmosheric screens should have a prefix, for example 'screen' and be 
        followed by an iterger. The default is 'screen'.
    
    Returns
    '''

    image=np.zeros((nPix,nPix))
    np.random.seed(0)
    amp= np.random.normal(0, SDIslandEffect, 6) #um
    #print(amp)
   
    #print("SD:\n")
    SD= statistics.stdev(amp)
    print(SD)
    name= Mod_segment(amp,1)
    hdu_tt= fits.open('Time_hist_TT_wind10ms_TelZ45M2_290sec_500Hz.fits')
    Data_tt = (hdu_tt [0] .data)*1e-09
    tip=(Data_tt[0]).astype(float)
    tilt=(Data_tt[1]).astype(float)
    SD_tip= statistics.stdev(tip)
    SD_tilt=statistics.stdev(tilt)
    Data_tip=tip*SD_WS/SD_tip
    Data_tilt=tilt*SD_WS/SD_tilt
    SD_tt_norm= statistics.stdev(Data_tip)
    SD_tt_norm= statistics.stdev(Data_tilt)
    
    for i in np.arange(1,nScreen):
        print('\r{:04}'.format(i),end = '')
        if mag==8 or mag==0:
            atmScreen = po.FITSOpticalElement( transmission=pupilFileName
                                           ,opd='C:/Users/esoria.DOMAINT/Desktop/ANDES/phase_screens_8M/part3-screen{:.0f}.fits'
                                           .format(i+100),opdunits='meters'
                                          ,planetype=po.poppy_core.PlaneType.pupil)
        if mag==12 or mag==14:
            atmScreen = po.FITSOpticalElement( transmission=pupilFileName
                                           ,opd='C:/Users/esoria.DOMAINT/Desktop/ANDES/phase_screens_12M/port3-screen{:.0f}.fits'
                                           .format(i+100),opdunits='meters'
                                          ,planetype=po.poppy_core.PlaneType.pupil)
        if mag==14:
            atmScreen.opd*=1.6
        if mag==0:
            atmScreen.opd*=0
        piston_segm=po.FITSOpticalElement(transmission='Tel-Pupil.fits', opd='C:/Users/esoria.DOMAINT/Desktop/ANDES/static/'+name.format(i)
                                          ,opdunits='nm',planetype=po.poppy_core.PlaneType.pupil)
        
        osys = po.OpticalSystem(npix = 400,oversample = 4)
        osys.add_pupil(atmScreen)
      
        if SD_WS!=0:
            wfe = po.ZernikeWFE(radius=19*u.m,
                                coefficients=[0, Data_tip[i+100], Data_tilt[i+100], 0e-6, 0, 0e-6, 0e-6],
                                aperture_stop=False)
            osys.add_pupil(wfe)
        if coro!=0:
            osys.add_image()
            osys.add_image(po.BandLimitedCoron(kind='cicular',  sigma=5.0))
        osys.add_pupil(piston_segm)
       
        amp= np.random.normal(0, SDSegmentJitter, 6) #um
        name2= Mod_segment(amp,2)
        piston_vibra=po.FITSOpticalElement(transmission='Tel-Pupil.fits', opd='C:/Users/esoria.DOMAINT/Desktop/ANDES/vibra/'+name2.format(i)
                                          ,opdunits='nm',planetype=po.poppy_core.PlaneType.pupil)
        osys.add_pupil(piston_vibra)
        if nPix==400:
            osys.add_detector(pixelscale=0.005, fov_arcsec=0.5)
        if nPix==800:
            osys.add_detector(pixelscale=0.005, fov_arcsec=1)
        # psf = osys.calc_psf(1.6e-6, display_intermediates=True)#for debugging
        psf = osys.calc_psf(wavelength)
        image += psf[0].data
        i=+1
    
    print('\nfinished')


    psf[0].data = image.copy()
    psf[0].header['MAGN'] = (mag ,'magnitud reference star')
    psf[0].header['SD_SP'] = (round(SD,3) ,'SD static piston (um)')
    psf[0].header['SD_DP'] = (SDSegmentJitter ,'SD dynamic piston (um)')
    psf[0].header['TT'] = (SD_tt_norm ,'SD TT(nm)')
    filename1 = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
    if mag==0 :
        name=filename1+'I={}-{:.0f}nm-DL.fits'.format(mag,wavelength*10**9)
    elif mag==8:
        name=filename1+'-I={}-{:.0f}nm.fits'.format(mag,wavelength*10**9)
    elif mag == 12: 
        name=filename1+'-I={}-{:.0f}nm.fits'.format(mag,wavelength*10**9)
    elif mag == 14:
        name=filename1+'-I={}-{:.0f}nm.fits'.format(mag,wavelength*10**9)

    psf.writeto(name,overwrite = True)
    print('sucessfull write')

    return name, SD


