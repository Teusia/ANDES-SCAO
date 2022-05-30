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
            ,atmWorseningFactor,seed,pupilMap,pupilFileName='Tel-Pupil.fits', atmFilePrefix='screen'
            ,outputFileSufix='I=8'):
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
    seed : float
        seed for the random noise generator
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
    -------
    imageFileName : str
        full name of the image generated at the end of this function.
    '''

    image=np.zeros((nPix,nPix))
    rng = np.random.default_rng(seed)
    amp= rng.normal(0, rmsIslandEffect, 6) #mm
    staticPist = Mod_segment(pupilMap,amp)
    
    for i in np.arange(1,nScreen):
        print('\r{:04}'.format(i),end = '')
       
        atmScreen = po.FITSOpticalElement( transmission=pupilFileName
                                       ,opd=atmFilePrefix+'{:.0f}.fits'
                                       .format(i),opdunits='meters'
                                      ,planetype=po.poppy_core.PlaneType.pupil)
        
        # tmp = atmScreen.opd.copy()
        
        atmScreen.opd*=atmWorseningFactor
        jit = rng.normal(0,rmsSegmentJitter,6)
        jitterPist = Mod_segment(pupilMap, jit)
        
        atmScreen.opd+=staticPist+jitterPist
        # piston_segm=po.FITSOpticalElement(transmission=pupilFileName, 
        #                                   opd='Segm_pupil.fits'
        #                                   ,opdunits='nm',
        #                                   planetype=po.poppy_core.PlaneType.pupil)
        
        osys = po.OpticalSystem(npix = 400,oversample = 4)
        osys.add_pupil(atmScreen)
       
       
       # wfe = po.ZernikeWFE(radius=15*u.m,
                        #    coefficients=[np.random.random()*1e-7, np.random.random()*1e-7, 1e-7, 0e-6, 0, 0e-6, 0e-6],
                          #  aperture_stop=True)
      
        #osys.add_pupil(wfe)
        
        osys.add_detector(pixelscale=0.005, fov_arcsec=0.5)
        # psf = osys.calc_psf(1.6e-6, display_intermediates=True)#for debugging
        psf = osys.calc_psf(wavelength)
        image += psf[0].data
    
    print('\nfinished')


    psf[0].data = image.copy()
    psf[0].header['ATMFACTO'] = (atmWorseningFactor ,'factor applied to the atmospheric screen')
    psf[0].header['RNGSEED'] = seed
    psf[0].header['RMSSTAT'] = rmsIslandEffect
    psf[0].header['RMSJIT'] = rmsSegmentJitter
    psf[0].header['MAGNITUD'] = outputFileSufix
    filename1 = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")+'_'+outputFileSufix+'_{:.0f}nm.fits'.format(wavelength*10**9)
    psf.writeto(filename1)
    
    print('sucessfull write')

    return filename1


if __name__ == '__main__':
    ppMap = buildPupil()
    
    outFile = SCAOSim(1600*10**-9, 400, 2000, 500*10**-9, 50*10**-9, 0, 1.0
                      , 12345, ppMap)




