# -*- coding: utf-8 -*-
"""
Created on Thu Mar  3 14:58:39 2022

@author: anne.cheffot
"""
from astropy.io import fits
import matplotlib.pyplot as plt

part = 3
startScreen = 200
for i in range(2000):
    with fits.open('Res-OPDseq-nm_PART-{:.0f}.fits'.format(part)) as hdul1:
        # hdul1.info()
        test = hdul1[0].data[i+startScreen,0,:,:]*10**-9
        
        wavelength = 850*10**-9
        telDia = 38.542
    
        
        plt.figure(1)
        plt.clf()
        plt.imshow(test)
        
        
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
        hdul.writeto('part{:.0f}-screen{:.0f}.fits'.format(part,startScreen+i)
                     ,overwrite = True)
    
