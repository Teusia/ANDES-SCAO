# -*- coding: utf-8 -*-
"""
Created on Fri Jan 28 11:45:35 2022

@author: anne.cheffot
"""
import numpy as np
import matplotlib.pyplot as plt
import poppy as po
from astropy.io import fits
from astropy import units as u

#if you want primary mirror segmentation use this. but not recomended
# seglist = [a for a in range((3*4**2)+(3*4)+1,(3*15**2)+(3*15)+1)]
# seglist.extend([a for a in range(723,736)])
# seglist.extend([a for a in range(739,752)])
# seglist.extend([a for a in range(755,768)])
# seglist.extend([a for a in range(771,784)])
# seglist.extend([a for a in range(787,800)])
# seglist.extend([a for a in range(803,816)])

# seglist.extend([a for a in range(821,831)])
# seglist.extend([a for a in range(838,848)])
# seglist.extend([a for a in range(855,865)])
# seglist.extend([a for a in range(872,882)])
# seglist.extend([a for a in range(889,899)])
# seglist.extend([a for a in range(906,916)])

# ap = po.MultiHexagonAperture(rings = 15, flattoflat = 1.239*u.m,gap=0.004*u.m,
#                               segmentlist = seglist,
#                               rotation = 0)
# pupilStop = po.CircularAperture(radius=38.542/2)
# sec = po.SecondaryObscuration(secondary_radius= 6.5*u.m, n_supports = 6,
#                               support_width = 0.54*u.m,rotation = 30)
# pupil = po.CompoundAnalyticOptic(opticslist=[
#     ap,pupilStop,
#     sec
#     ],name = 'ELT pupil')

image = np.zeros((400,400))
wavelength = 2.2*10**-6
startScreen = 200
# worseningFactor = 0.0#diffraction limited
worseningFactor = 1.0 #I = 8
# worseningFactor = 1.6 #I = 11.5
# worseningFactor = 2.16 #I = 15
# worseningFactor = 2.57 #I = 14


for i in range(2000):
    print('\r{:04}'.format(i),end = '')
    
    
    atmScreen = po.FITSOpticalElement( transmission='Tel-Pupil.fits'
                                       ,opd='part3-screen{:.0f}.fits'
                                       .format(startScreen+i)
                                      ,planetype=po.poppy_core.PlaneType.pupil)
    
    tmp = atmScreen.opd.copy()
    
    atmScreen.opd*=worseningFactor
    
    osys = po.OpticalSystem(npix = 400,oversample = 4)
    osys.add_pupil(atmScreen)
    osys.add_detector(pixelscale=0.005, fov_arcsec=0.5)
    # psf = osys.calc_psf(1.6e-6, display_intermediates=True)#for debugging
    psf = osys.calc_psf(wavelength)
    image += psf[0].data

print('\nfinished')


psf[0].data = image.copy()
psf[0].header['ATMFACTO'] = (worseningFactor,'factor applied to the atmospheric screen')
if worseningFactor==0.0 :
    psf.writeto('part3-2000screen-I=8-{:.0f}nm-diffractionLimited.fits'.format(wavelength*10**9),overwrite = True)
elif worseningFactor==1.0 :
    psf.writeto('part3-2000screen-I=8-{:.0f}nm.fits'.format(wavelength*10**9),overwrite = True)
elif worseningFactor == 1.6: 
    psf.writeto('part3-2000screen-I=11_5-{:.0f}nm-update.fits'.format(wavelength*10**9),overwrite = True)
elif worseningFactor == 2.16:
    psf.writeto('part3-2000screen-I=15-{:.0f}nm-core.fits'.format(wavelength*10**9),overwrite = True)
elif worseningFactor == 2.57:
    psf.writeto('part3-2000screen-I=14-{:.0f}nm-core.fits'.format(wavelength*10**9),overwrite = True)
print('sucessfull write')


