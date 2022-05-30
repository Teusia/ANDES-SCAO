"""
Created on Tue May 17 16:30:42 2022

@author: EstherSoria
"""
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
#inputs must be an array

def buildPupil(pupilFileName='Tel-Pupil.fits'):
    '''
    build the pupil map to allow adressing single segments
    
    Parameters
    ----------
    pupilFileName : string
        name of the .fits file that contain the pupil shape

    Returns
    -------
    pupMap : array of int
        2D map of the pupil. each pixels that belong to a segment is numbered 
        after the segment number

    '''
    
    hdu = fits.open(pupilFileName)
    Data1 = hdu[0].data # get data
    
    #generate a meshgrid with the pixel position to evaluate the position with respect its angle
    x = np.linspace(-1,1,Data1.shape[0])
    X,Y = np.meshgrid(x,x)
    theta = np.arctan2(Y,X)#I am interested only in the angle each pixel 
    limites = np.arange(-5*np.pi/6,np.pi,2*np.pi/6)
    
    pupMap = Data1.copy()
    tmppupilMask = Data1.copy()
    #generate a pupil map with the number of the segment in which the pixel is placed
    for s in range(5):
         mask = np.logical_and(theta[:,:]>limites[s] , theta[:,:]<limites[s+1]) #codition     
         mask = mask*tmppupilMask
         pupMap += mask*(s+1)       
    
    return pupMap
    

def Mod_segment(pMap,amp):
    '''
    create piston errors of each segments

    Parameters
    ----------
    pMap : array of int
        pupil map as created by buildPupil
    amp : list of float
        amplitude of the piston for each segments.

    Returns
    -------
    newPupil : array of float
        phase error of the pupil.

    '''
    
    newPupil=np.zeros(pMap.shape)
    amp= np.asarray(amp) #convert into nm
    #Give to each segment the value of the amplitud
    
    #applying the value of the amplitudes
    for s in range(6):
         newPupil[pMap==s+1] = amp[s] 
         
    return newPupil
