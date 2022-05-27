"""
Created on Tue May 17 16:30:42 2022

@author: EstherSoria
"""
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
#inputs must be an array
def Mod_segment(amp,types):
    hdu = fits.open('Tel-Pupil.fits')
    Data1 = hdu [0] .data # get data
    cx=Data1//2 #center of the pupil x coordinate
    cy=Data1//2 #center of the pupil y coordinate
   
    #Give to each segment the value of the amplitud
    Newpupil=Data1.astype('float64')
    #generate a meshgrid with the pixel position to evaluate the position with respect its angle
    x = np.linspace(-1,1,Data1.shape[0])
    X,Y = np.meshgrid(x,x)
    theta = np.arctan2(Y,X)#I am interested only in the angle each pixel 
    limites = np.arange(-5*np.pi/6,np.pi,2*np.pi/6)
    
    pupilMap = Data1.copy()
    tmppupilMask = Data1.copy()
    #generate a pupil map with the number of the segment in which the pixel is placed
    for s in range(5):
         mask = np.logical_and(theta[:,:]>limites[s] , theta[:,:]<limites[s+1]) #codition     
         mask = mask*tmppupilMask
         pupilMap += mask*(s+1)       
        
    #applying the value of the amplitudes
    for s in range(6):
         Newpupil[pupilMap==s+1] = amp[s] 
         
    #Saving in poppy format
    part = 3
    wavelength = 850*10**-9
    telDia = 38.542
    hdu = fits.PrimaryHDU(Newpupil)
    hdu.header['WAVELEN'] = wavelength
    hdu.header['DIAM'] = telDia
    hdu.header['DIFFLMT'] = wavelength/telDia
    hdu.header['OVERSAMP'] = 1
    hdu.header['DET_SAMP'] = 1
    hdu.header['PUPLSCAL'] = telDia/Data1.shape[0] 
    hdu.header['PIXUNIT'] = 'meters'
    hdu.header['BUNIT'] = 'meters'
    hdul = fits.HDUList(hdu)
    if types==1:
        hdul.writeto('Segm_pupil.fits'.format(part),overwrite = True)
        plt.imshow(Newpupil)
        plt.colorbar(label="$n$m") 
        plt.show() 
    if types==2:
         hdul.writeto('Vibra_pupil.fits'.format(part),overwrite = True)      
    return 
