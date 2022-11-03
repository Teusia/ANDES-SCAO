# -*- coding: utf-8 -*-
"""
Created on Fri Jan 28 11:45:35 2022

@author: anne.cheffot
"""
import numpy as np
import matplotlib.pyplot as plt
import poppy as po
from astropy.io import fits

with fits.open('Res-OPDseq-nm_PART-3.fits') as hdul:
    hdul.info()
    test = hdul[0].data[0,0,:,:]

seglist = [a for a in range((3*4**2)+(3*4)+1,(3*15**2)+(3*15)+1)]
seglist.extend([a for a in range(723,736)])
seglist.extend([a for a in range(739,752)])
seglist.extend([a for a in range(755,768)])
seglist.extend([a for a in range(771,784)])
seglist.extend([a for a in range(787,800)])
seglist.extend([a for a in range(803,816)])

seglist.extend([a for a in range(821,831)])
seglist.extend([a for a in range(838,848)])
seglist.extend([a for a in range(855,865)])
seglist.extend([a for a in range(872,882)])
seglist.extend([a for a in range(889,899)])
seglist.extend([a for a in range(906,916)])

ap = po.MultiHexagonAperture(rings = 15, flattoflat = 1.239,gap=0.004,
                              segmentlist = seglist,
                              rotation = 0)
pupilStop = po.CircularAperture(radius=38.542/2)
sec = po.SecondaryObscuration(secondary_radius= 6.5, n_supports = 6,
                              support_width = 0.54,rotation = 30)
pupil = po.CompoundAnalyticOptic(opticslist=[
    ap,
    pupilStop,
    sec
    ],name = 'ELT pupil')
# pupil.display(npix = 400)
# plt.figure(6)
# plt.clf()
# fig,axs2 = 
# plt.imshow(pupil.transmission)
#%%
image = np.zeros((400,400))
wavelength = 1.6*10**-6
spart = 3
startScreen = 200
# worseningFactor = 1.0 #I = 8
worseningFactor = 1.6 #I = 11.5
# worseningFactor = 2.16 #I = 15
# worseningFactor = 2.57#I=14

wavelength = 850*10**-9
telDia = 38.542
fourierMax = 1/(2*telDia/400)
x = np.linspace(-fourierMax, fourierMax,400)
X,Y = np.meshgrid(x,x)
R = np.sqrt(X**2+Y**2)

for i in range(1):
    print('\r{:04}'.format(i),end = '')
    
    
    atmScreen = po.FITSOpticalElement( transmission='Tel-pupil.fits'
                                      ,opd='part{:.0f}-screen{:.0f}.fits'
                                      .format(spart,startScreen+i)
                                      ,planetype=po.poppy_core.PlaneType.pupil)
    tmp = atmScreen.opd.copy()
    plt.figure(10)
    plt.clf()
    fig,axs = plt.subplots(ncols = 3, num = 10)
    axs[0].imshow(tmp,vmin=-10**-7,vmax = 10**-7)
    plt.figure(11)
    plt.clf()
    fig1,axs1 = plt.subplots(ncols = 2,num = 11)
    
    fourierSpace = np.fft.fftshift(np.fft.fft2(tmp))*tmp.shape[0]
    axs1[0].imshow(np.abs(fourierSpace)**2)
    fourierSpace[R<1]*=worseningFactor
    axs1[1].imshow(np.abs(fourierSpace)**2)
    
    pupilSpace = np.abs(np.fft.ifft2(fourierSpace)/tmp.shape[0])
    pupilSpace *= atmScreen.amplitude
    pupilSpace[atmScreen.amplitude.astype(bool)] -= np.mean(pupilSpace[atmScreen.amplitude.astype(bool)])
    axs[1].imshow((pupilSpace),vmin=-10**-7,vmax = 10**-7)
    axs[2].imshow(tmp-(pupilSpace),vmin=-10**-7,vmax = 10**-7)
    atmScreen.opd = pupilSpace.copy()
    
    # atmScreen.opd*=worseningFactor
    
    osys = po.OpticalSystem(npix = 400,oversample = 4)
    # osys = po.OpticalSystem(oversample = 4)
    osys.add_pupil(pupil)
    osys.add_pupil(atmScreen)
    osys.add_detector(pixelscale=0.005, fov_arcsec=0.5)
    # psf = osys.calc_psf(1.6e-6, display_intermediates=True)
    psf = osys.calc_psf(wavelength)
    image += psf[0].data

    plt.figure(2)
    plt.clf()
    plt.imshow(np.log(image))
    # plt.pause(0.1)

# plt.plot(4.75*np.cos(np.linspace(-1,1,51)*np.pi)+199.5,
#          4.75*np.sin(np.linspace(-1,1,51)*np.pi)+199.5)
# firstExtinction = psf[0].header['DIFFLMT']*1.22/psf[0].header['PIXELSCL']
# plt.plot(firstExtinction*np.cos(np.linspace(-1,1,51)*np.pi)+199.5,
#          firstExtinction*np.sin(np.linspace(-1,1,51)*np.pi)+199.5)

psf[0].data = image.copy()
psf[0].header['ATMFACTO'] = (worseningFactor,'factor applied to the atmospheric screen')
# psf.writeto('part3-2000screen-I=15-{:.0f}nm.fits'.format(wavelength*10**9))
# hdu = fits.PrimaryHDU(image)
# hdul = fits.HDUList(hdu)
# hdul.writeto('part3-2000screen-1photeach.fits')
# plt.imshow(np.log(psf[0].data))

# plt.imshow((psf[0].data))
# po.display_psf(psf, title = 'elt psf')

