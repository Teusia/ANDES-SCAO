"""
Created on Fri Jan 28 11:45:35 2022

@author: anne.cheffot
@author: EstherSoria

"""
import numpy as np
import poppy as po
import datetime
from pupil_ANDES import *

def SCAOSim(wavelength,nPix,nScreen,rmsIslandEffect,rmsSegmentJitter,rmsWindshake
            ,atmWorseningFactor,pupilFileName='Tel-Pupil.fits', atmFilePrefix='screen'
            ,outputFileSufix='I=8-2200nm.fits'):
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

    image=np.zeros((nPix,nPix))
    amp= np.random.normal(0, rmsIslandEffect, 6) #um
    Mod_segment(amp,1)
    
    for i in np.arange(1,nScreen):
        print('\r{:04}'.format(i),end = '')
       
        atmScreen = po.FITSOpticalElement( transmission=pupilFileName
                                       ,opd='part3-screen{:.0f}.fits'
                                       .format(i),opdunits='meters'
                                      ,planetype=po.poppy_core.PlaneType.pupil)
        
        tmp = atmScreen.opd.copy()
        
        atmScreen.opd*=atmWorseningFactor
        piston_segm=po.FITSOpticalElement(transmission='Tel-Pupil.fits', opd='Segm_pupil.fits'.format(i)
                                          ,opdunits='nm',planetype=po.poppy_core.PlaneType.pupil)
        
        osys = po.OpticalSystem(npix = 400,oversample = 4)
        osys.add_pupil(atmScreen)
       
        amp= np.random.normal(0, rmsSegmentJitter, 6) #um
        Mod_segment(amp,2)
       
       # wfe = po.ZernikeWFE(radius=15*u.m,
                        #    coefficients=[np.random.random()*1e-7, np.random.random()*1e-7, 1e-7, 0e-6, 0, 0e-6, 0e-6],
                          #  aperture_stop=True)
      
        #osys.add_pupil(wfe)
        
        osys.add_pupil(piston_segm)
        piston_vibra=po.FITSOpticalElement(transmission='Tel-Pupil.fits', opd='Vibra_pupil.fits'.format(i)
                                          ,opdunits='nm',planetype=po.poppy_core.PlaneType.pupil)
        osys.add_pupil(piston_vibra)
        osys.add_detector(pixelscale=0.005, fov_arcsec=0.5)
        # psf = osys.calc_psf(1.6e-6, display_intermediates=True)#for debugging
        psf = osys.calc_psf(wavelength)
        image += psf[0].data
    
    print('\nfinished')


    psf[0].data = image.copy()
    psf[0].header['ATMFACTO'] = (atmWorseningFactor ,'factor applied to the atmospheric screen')
    filename1 = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
    if atmWorseningFactor==0.0 :
        psf.writeto(filename1+'I=8-{:.0f}nm-DL.fits'.format(wavelength*10**9),overwrite = True)
    elif atmWorseningFactor==1.0 :
        psf.writeto(filename1+'-I=8-{:.0f}nm.fits'.format(wavelength*10**9),overwrite = True)
    elif atmWorseningFactor == 1.6: 
        psf.writeto(filename1+'-I=11_5-{:.0f}nm.fits'.format(wavelength*10**9),overwrite = True)
    elif atmWorseningFactor == 2.16:
        psf.writeto(filename1+'-I=15-{:.0f}nm.fits'.format(wavelength*10**9),overwrite = True)
    elif atmWorseningFactor == 2.57:
        psf.writeto(filename1+'-I=14-{:.0f}nm.fits'.format(wavelength*10**9),overwrite = True)
    print('sucessfull write')

    return
