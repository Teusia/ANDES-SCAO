"""
Created on Tue May 17 16:30:42 2022

@author: EstherSoria
"""
from astropy.io import fits
import matplotlib.pyplot as plt
import cv2
import numpy as np
import numpy as np
#inputs must be an array
def Mod_segment(amp,segm,types):
    hdu = fits.open('Tel-Pupil.fits')
    Data1 = hdu [0] .data # get data
    cx=200 #center of the pupil x coordinate
    cy=200 #center of the pupil y coordinate
    ampli=np.ones((6,1))
    #create the vectors to allocate the idex of each segment:
    """
    for i in np.arange(1,7,1):
        u=exec("gx%s = []" % i)
        v=exec ("gy%s = []" % i)
        globals()[u] =[]
        globals()[v] =[]
    """   
    gx1=[]; gx2=[]; gx3=[]; gx4=[]; gx5=[]; gx6=[]
    gy1=[]; gy2=[]; gy3=[]; gy4=[]; gy5=[]; gy6=[]

    #evaluate the postition of each ON-pixel with respect to its angle:
    for i in np.arange(0,np.size(Data1,0)):       
        for j in np.arange(0,np.size(Data1,1)):
            pixel=Data1[i,j]
            if pixel==1:       
                theta=np.arctan2(j-cx,i-cy)
                if theta>2.1 and theta<3.15:
                    gx1.append(i)
                    gy1.append(j)
                if theta>1.05 and theta<2.1:
                    gx2.append(i)
                    gy2.append(j)
                if theta>0 and theta<1.05:
                    gx3.append(i)
                    gy3.append(j)
                if theta>-1.05 and theta<0:
                    gx4.append(i)
                    gy4.append(j)
                if theta>-2.1 and theta<-1.05:
                    gx5.append(i)
                    gy5.append(j)
                if theta>-3.15 and theta<-2.1:
                    gx6.append(i)
                    gy6.append(j)
            else:
                pass
   
    #Give to each segment the value of the amplitud
    Newpupil=Data1.astype('float64')
    u = np.size(segm)
    for i in range(u):
        ampli[segm[i]-1]=amp[segm[i]-1]
    Newpupil[gx1,gy1]=ampli[0]*10**3
    Newpupil[gx2,gy2]=ampli[1]*10**3
    Newpupil[gx3,gy3]=ampli[2]*10**3
    Newpupil[gx4,gy4]=ampli[3]*10**3
    Newpupil[gx5,gy5]=ampli[4]*10**3
    Newpupil[gx6,gy6]=ampli[5]*10**3
    
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
    if types==1: #simulating piston error between segments
        hdul.writeto('Segm_pupil.fits'.format(part),overwrite = True)
        plt.imshow(Newpupil)
        plt.colorbar(label="$n$m") 
        plt.show() 
    if types==2: #simulating vibrations
         hdul.writeto('Vibra_pupil.fits'.format(part),overwrite = True)      
    return ampli
