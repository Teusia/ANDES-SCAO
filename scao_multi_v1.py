# -*- coding: utf-8 -*-
"""
Created on Mon Jun  6 13:02:21 2022

@author: esoria
"""

"""
Created on Fri Jan 28 11:45:35 2022

@author: anne.cheffot
@author: EstherSoria
"""
import numpy as np
import poppy as po
import datetime
import statistics
from pupil_ANDES import *
from sklearn.metrics import mean_squared_error
import math
import time
from astropy import units as u

def SCAOSim_mul(mag,wavelength,value):
    '''
    Parameters
    ----------
    wavelength : float
        wavelength of the final image in meters
    nPix : int
        number of pixels across one image
    nScreen : int
        number of time light is propagated through an aberated screen.
    rmsIslandEffect : float
        root mean square error of the island effect in meters
    rmsSegmentJitter : float
        root mean square of the segment jitter in meters
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
    outputFileSufix : string, optional
        The output file should be date_time followed by a suffix.
        The default is 'I=8-2200nm.fits'.
    Returns
    '''
    nPix=400
    nScreen=200
    SDSegmentJitter=0
    rmsWindshake=0
    pupilFileName='Tel-Pupil.fits'
    atmFilePrefix='screen'
    outputFileSufix='I=8-2200nm.fits'
    image=np.zeros((nPix,nPix))
    #amp= np.random.normal(0, SDIslandEffect, 6) #um
    #np.random.shuffle(amp)
    amp=np.asarray([-0.09718888,  0.0276449 ,  0.01844131,  0.01865384,  0.11139951,
        0.02556878])
    amp=amp*value
    print(amp)
    '''
    y_predicted = [0,0,0,0,0,0]
    MSE = mean_squared_error(amp, y_predicted)
    RMSE = math.sqrt(MSE)
    print("Root Mean Square Error:\n")
    print(RMSE)
    '''
    SD= statistics.stdev(amp)
    name=Mod_segment(amp,1)
    
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
        if rmsWindshake!=0:
            wfe = po.ZernikeWFE(radius=19*u.m,
                                coefficients=[0, np.random.random()*rmsWindshake, np.random.random()*rmsWindshake, 0e-6, 0, 0e-6, 0e-6],
                                aperture_stop=False)
            osys.add_pupil(wfe)
        
        osys.add_pupil(piston_segm)
        amp= np.random.normal(0, SDSegmentJitter, 6) #um
        #Mod_segment(amp,2)
        #piston_vibra=po.FITSOpticalElement(transmission='Tel-Pupil.fits', opd='Vibra_pupil.fits'.format(i)
                                       #   ,opdunits='nm',planetype=po.poppy_core.PlaneType.pupil)
        #osys.add_pupil(piston_vibra)
        osys.add_detector(pixelscale=0.005, fov_arcsec=0.5)
        # psf = osys.calc_psf(1.6e-6, display_intermediates=True)#for debugging
        psf = osys.calc_psf(wavelength)
        image += psf[0].data
        i=+1
    
    print('\nfinished')


    psf[0].data = image.copy()
    psf[0].header['MAGN'] = (mag ,'magnitud reference star')
    psf[0].header['SD_SP'] = (round(SD,3) ,'SD static piston (um)')
    psf[0].header['SD_DP'] = (SDSegmentJitter ,'SD dynamic piston (um)')
    psf[0].header['TT'] = (rmsWindshake ,'TT error (m)')
    filename1 = datetime.datetime.now().strftime("%Y%m%d-%H%M%S%f")
    if mag==0 :
        name=filename1+'I=0-{:.0f}nm-DL.fits'.format(wavelength*10**9)
    elif mag==8:
        name=filename1+'-I=8-{:.0f}nm.fits'.format(wavelength*10**9)
    elif mag == 12: 
        name=filename1+'-I=12-{:.0f}nm.fits'.format(wavelength*10**9)
    elif mag == 14:
        name=filename1+'-I=14-{:.0f}nm.fits'.format(wavelength*10**9)
   
    psf.writeto(name,overwrite = True)
    print('sucessfull write')

    return name, SD


