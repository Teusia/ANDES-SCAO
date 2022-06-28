# -*- coding: utf-8 -*-
"""
Created on Tue Mar  8 17:49:34 2022

@author: anne.cheffot
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import scipy.optimize as opt
from matplotlib.lines import Line2D
import os
from cog import cog


plt.rcParams.update({'font.size': 15})

def computeRadialProfile(image, centerInPxY, centerInPxX):
    yCoord, xCoord= np.indices(image.shape)
    yCoord= (yCoord - centerInPxY)
    xCoord= (xCoord - centerInPxX)
    rCoord=np.sqrt(xCoord**2 + yCoord**2)
    indexR= np.argsort(rCoord.flat)
    radialDistancesSorted= rCoord.flat[indexR]
    imageValuesSortedByRadialDistance= image.flat[indexR]
    integerPartOfRadialDistances= radialDistancesSorted.astype(int)
    deltaRadialDistance= integerPartOfRadialDistances[1:] - \
        integerPartOfRadialDistances[:-1]
    radialDistanceChanges= np.where(deltaRadialDistance)[0]
    nPxInBinZero= radialDistanceChanges[0]+ 1
    nPxInRadialBin= radialDistanceChanges[1:] - \
        radialDistanceChanges[:-1]
    imageRadialCumSum= np.cumsum(imageValuesSortedByRadialDistance,
                                 dtype=np.float64)
    imageSumInBinZero= imageRadialCumSum[radialDistanceChanges[0]]
    imageSumInRadialBin= \
        imageRadialCumSum[radialDistanceChanges[1:]] - \
        imageRadialCumSum[radialDistanceChanges[:-1]]
    profileInZero= imageSumInBinZero / nPxInBinZero
    profileFromOne= imageSumInRadialBin / nPxInRadialBin
    profile= np.hstack([profileInZero, profileFromOne])

    distanceRadialCumSum= np.cumsum(radialDistancesSorted)
    distanceSumInBinZero= distanceRadialCumSum[radialDistanceChanges[0]]
    distanceSumInRadialBin= \
        distanceRadialCumSum[radialDistanceChanges[1:]] - \
        distanceRadialCumSum[radialDistanceChanges[:-1]]
    distanceInZero= distanceSumInBinZero / nPxInBinZero
    distanceFromOne= distanceSumInRadialBin / nPxInRadialBin
    radialDistance= np.hstack([distanceInZero, distanceFromOne])
    return profile, radialDistance

def twoD_Gaussian(coordTuple, amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
    (x, y) = coordTuple
    xo = float(xo)
    yo = float(yo)    
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) 
                            + c*((y-yo)**2)))
    return g.ravel()

def adjustCenter(img,dbg = False):
    xsize = img.shape[0]
    xx,yy = np.meshgrid(np.arange(xsize),np.arange(xsize))
    
    popt, pcov = opt.curve_fit(twoD_Gaussian, (xx.ravel(), yy.ravel()), 
                               img.ravel(), p0=(np.max(img),xsize/2,
                                               xsize/2,1,1,0,0))
    if dbg:
        plt.figure('debug adjustCenter')
        plt.clf()
        plt.imshow(np.log(img))
        fitted = twoD_Gaussian((xx.ravel(),yy.ravel()),*popt)
        plt.contour(fitted.reshape(xsize,xsize))
    return popt[1], popt[2]

def binnedContrast(psfData,binSize,pix2arcsec):
    '''
    this function takes the path to the file of the psf and file name, and the 
    bin size in arcsec
    
    parameters
    ----------
    psf: str
        path and file name of the PSF.
    binSize : float
        size of the bin in mas
    
    output
    ------
        binCoord: array of float
            position of the bin center including the center fo the PSF
        binAvr: average of each bin distance
        binMin: minimum of each bin ditance
        binMax: maximum of each bin ditance
        binStd: standard deviation of each bin
    
    '''
    #canculate the field of view
    fov = psfData.shape[0]*pix2arcsec
    #calculate the number of bin and their limites
    nbBin = np.round(fov/2/binSize).astype(int)
    binEdges = np.arange(binSize/2,fov/2,binSize)
    
    #calculate the reference
    mask = R<=binSize/2
    ref = np.sum(test*mask)
    
    #initialise variables to fill
    mmin = np.zeros(nbBin)
    mmax= np.zeros(nbBin)
    maverage = np.zeros(nbBin)
    mstd = np.zeros(nbBin)
    mmin[0] = mmax[0] = maverage[0] = ref
    
    for rd in range(1,nbBin):#for each bin
        alpha = np.arctan(1/(2*rd))#find the displacement
        nbPos = np.round(2*np.pi/alpha).astype(int)#the number of displacement
        curring = []
        for ap in range(nbPos):
            #build the mask to select the wanted pixel
            tmpX,tmpY = np.meshgrid(x-(rd*binSize*np.cos(alpha*ap)),
                                    x-(rd*binSize*np.sin(alpha*ap)))
            tmpR = np.sqrt(tmpX**2+tmpY**2)
            mask=tmpR<=binSize/2
            
            #accumulate the pixels sum
            curring.append(np.sum(test*mask))
        # plt.figure(10+rd)
        # plt.clf()
        # plt.imshow(np.log(test),cmap = 'hot',extent = [-200*pix2arcsec,200*pix2arcsec
        #                                        ,-200*pix2arcsec,200*pix2arcsec])
        # plt.scatter(rd*binSize*np.cos(alpha*np.arange(nbPos))
        #             ,rd*binSize*np.sin(alpha*np.arange(nbPos)),
        #             c=np.log(curring),cmap='hot')
        # plt.title('radial distance {:.0f} mas'.format(radDist*rd))

        # plt.figure(10+rd)
        # plt.clf()
        # plt.hist(curring/ref)
        # plt.title('radial distance {:.0f} mas'.format(radDist*rd))
        # plt.xlabel('contrast values')
        # plt.ylabel('number of samples')
        mmin[rd] = np.min(curring)
        mmax[rd] = np.max(curring)
        maverage[rd] = np.mean(curring)
        mstd[rd] = np.std(curring)
    
    mmin *=1/ref
    mmax *=1/ref
    maverage *=1/ref
    mstd *= 1/ref
    xvar = np.arange(0,binEdges[-1],binSize)
    
    return xvar,mmin,mmax,maverage,mstd

def plotBinnedContrast(binDat,wavelength,starType,saveFig = False,
                       figName = None,xlimp = [],ylimp = []):
    '''
    plot the binned contrast for you

    Parameters
    ----------
    binDat : tuple of numpy array
        preferably directly the output of binnedContrast function
    saveFig : bool (optional)
        set to True if you want the figure saved automatically
    figName : str (optionnal)
        Set the name for the figure saved. This is not optionnal if you have 
        set saveFig to True

    Returns
    -------
    None.

    '''
    if saveFig and (figName is None):
        print('give a name to the figure you want to save')
        return 0
    fig = plt.figure(figName)
    plt.clf()
    plt.semilogy(binDat[0]*10**3,binDat[1],drawstyle = 'steps-mid')
    plt.plot(binDat[0]*10**3,binDat[2],drawstyle = 'steps-mid')
    plt.errorbar(binDat[0]*10**3, binDat[3], yerr = binDat[4],capsize=5
                        ,drawstyle = 'steps-mid')
    
    plt.xlabel('distance from PSF center [mas]')
    plt.ylabel('relative flux')
    
    plt.title('wl = {:.0f} nm, {}'.format(wavelength,starType))
    custom_lines = [Line2D([0],[0],color='C0'),
                    Line2D([0],[0],color='C1'),
                    Line2D([0],[0],color='C2'),
                    ]
    plt.legend(custom_lines,['min','max','average+std',
                             ],
               handletextpad=0.2,
               labelspacing = 0.2)
    if ylimp:
        plt.ylim(ylimp)
    if xlimp:
        plt.xlabel(xlimp)
    plt.grid()
    plt.tight_layout()
    
    if saveFig:
        plt.savefig(figName)
    

#%%
path2fits = '20220530-162324/'
listFits = os.listdir(path2fits)
listFits = [a for a in listFits if '.fits' in a]

baseline = {1000:[[[25,85],[35,95]],[[10**-2,1.6*10**-3],[10**-2,1.6*10**-3]]],
            1600:[[[25,85],[35,95]],[[10**-2,1.0*10**-3],[10**-2,1.0*10**-3]]],
            2200:[[[25,85],[35,95]],[[10**-6,1.0*10**-6],[10**-6,1.0*10**-6]]]}
goal = {1000: [[[25,85],[35,95]],[[10**-3,2.0*10**-4],[10**-3,2.0*10**-4]]],
        1600: [[[25,85],[35,95]],[[10**-3,1.0*10**-4],[10**-3,1.0*10**-4]]],
        2200:[[[25,85],[35,95]],[[10**-6,1.0*10**-6],[10**-6,1.0*10**-6]]]}
ylimplot = [9*10**-6,1.5]



plt.close('all')

radialStop = 0.25 #arcsec
binSize = 0.010#arcsec

magnitude = ['diffraction limit',8,15,14]
listCases = []

for lf in listFits:

    with fits.open(path2fits+lf) as hdul1:
        # hdul1.info()
        test = hdul1[0].data
        wl = int(hdul1[0].header['WAVELEN']*10**9)
        mag = hdul1[0].header['MAGNITUD']
        listCases.append([wl,mag])
        pix2arcsec = hdul1[0].header['PIXELSCL']
        hdul1.close()
    
    #normalise the image
    test *= 1/np.max(test)

    imgSize = test.shape[0]

    nbBin = np.floor(radialStop/binSize).astype(int)
    # binEdges = np.arange(binSize/2,radialStop,binSize)
    
    x = np.arange(-imgSize/2,imgSize/2)+0.5
    
    x*=pix2arcsec
    X,Y = np.meshgrid(x,x)
    R = np.sqrt(X**2+Y**2)
    
    binData = binnedContrast(test, binSize,pix2arcsec)
    plotBinnedContrast(binData,wl,mag,figName = lf)
    
    #adjust the psf center
    # centering = adjustCenter(test)
    centering = cog(test,0.05*np.max(test))
    
    #Calculate the psf profiles
    profile,radDist = computeRadialProfile(test, *centering)
    #Normalise the psf profile 
    nProfile = profile/np.max(profile)
    
    plt.figure(wl)
    mas = radDist*pix2arcsec*10**3
    plt.semilogy(mas,nProfile,label = mag)
    
    
    
#plots prettyfication
#psf profile
wluni = np.unique([a[0] for a in listCases])

for i in wluni:
    plt.figure(i)
    plt.plot(baseline[i][0],baseline[i][1],'r',linewidth = 3)
    plt.plot(goal[i][0],goal[i][1],'m',linewidth = 3)
    plt.grid()
    # custom_lines = [Line2D([0],[0],color='C0'),
    #                 Line2D([0],[0],color='C1'),
    #                 Line2D([0],[0],color = 'r'),
    #                 Line2D([0],[0],color = 'm')
    #                 ]
    plt.ylim(ylimplot)
    plt.legend(
        # custom_lines,['diffraction limit','I=8'
        #                        ,'Baseline','goal'
        #                      ],
            handletextpad=0.2,
                   labelspacing = 0.2)
    
    plt.xlim([0,100])
    plt.grid()
    plt.xlabel('distance from PSF center [mas]')
    plt.ylabel('relative intensity')
    
    plt.grid()
    plt.title('wavelength {:.0f} nm'.format(i))
    plt.tight_layout()

#binned psf profiles details of everything
# for lf in range(len(listFits)):
    
#     plt.figure(listFits[lf])
#     plt.plot(baseline[listCases[lf][0]][0],baseline[listCases[lf][0]][1],'r',linewidth = 3)
#     plt.plot(goal[listCases[lf][0]][0],goal[listCases[lf][0]][1],'m',linewidth = 3)
#     plt.xlabel('distance from PSF center [mas]')
#     plt.ylabel('relative flux')
    
#     plt.title('wl = {:.0f} nm, {}'.format(listCases[lf][0],listCases[lf][1]))
#     custom_lines = [Line2D([0],[0],color='C0'),
#                     Line2D([0],[0],color='C1'),
#                     Line2D([0],[0],color='C2'),
#                     ]
#     plt.legend(custom_lines,['min','max','average+std',
#                              ],
#                handletextpad=0.2,
#                labelspacing = 0.2)
#     plt.ylim(ylimplot)
#     plt.grid()
    # plt.tight_layout()

# plt.semilogy([90,90], [5*10**-5,1],'k--')
# plt.semilogy([30,30], [5*10**-5,1],'k--')