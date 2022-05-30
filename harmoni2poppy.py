# -*- coding: utf-8 -*-
"""
Created on Thu Mar  3 14:58:39 2022

@author: anne.cheffot
"""
from astropy.io import fits
import matplotlib.pyplot as plt

part = 3
startScreen = 0

for i in range(1000):
    # with fits.open('Res-OPDseq-nm_PART-{:.0f}.fits'.format(part)) as hdul1:
    with fits.open('resWF_13.fits'.format(part)) as hdul1:
        # hdul1.info()
        # test = hdul1[0].data[i+startScreen,0,:,:]*10**-9
        test = hdul1[0].data[i+startScreen,:,:]
        
        wavelength = 850*10**-9
        telDia = 38.542
        
        
        hdu = fits.PrimaryHDU(test)
        hdu.header['WAVELEN'] = wavelength
        hdu.header['DIAM'] = telDia
        hdu.header['DIFFLMT'] = wavelength/telDia
        hdu.header['OVERSAMP'] = 1
        hdu.header['DET_SAMP'] = 1
        hdu.header['PIXELSCL'] = telDia/test.shape[0] 
        hdu.header['PIXUNIT'] = 'meters'
        hdu.header['BUNIT'] = 'meters'
        hdul = fits.HDUList(hdu)
        hdul.writeto('screen_I=12/screenI=12{:.0f}.fits'.format(i+1000)
                     ,overwrite = True)
        
        hdul1.close()
    
