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
import json

plt.rcParams.update({'font.size': 15})

def computeRadialProfile(image, centerInPxY, centerInPxX):
    yCoord, xCoord= np.indices(image.shape)
    yCoord= (yCoord - centerInPxY)
    xCoord= (xCoord - centerInPxX)
    rCoord=np.sqrt(xCoord**2 + yCoord**2)
    indexR= np.argsort(rCoord.flat)
    radialDistancesSorted= rCoord.flat[indexR]
    imageValuesSortedByRadialDistance= image.flat[indexR]
    integerPartOfRadialDistances= radialDistancesSorted.astype(np.int)
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

waveNb = 0
wl = [1000,1600,2200]
fileList = [[
            # 'part3-2000screen-I=8-1000nm.fits',
            'part3-2000screen-I=8-1000nm-diffractionLimited.fits',
            'part3-2000screen-I=8-1000nm-update.fits',
            # 'part3-2000screen-I=11_5-1000nm.fits',
            # 'part3-2000screen-I=11_5-1000nm-update.fits'
            ],
            [
            # 'part3-2000screen-I=8-1600nm.fits',
            'part3-2000screen-I=8-1600nm-diffractionLimited.fits',
            'part3-2000screen-I=8-1600nm-update.fits',
            # 'part3-2000screen-I=11_5-1600nm.fits',
            # 'part3-2000screen-I=11_5-1600nm-update.fits'
            ],
             [
             # 'part3-2000screen-I=8-2200nm.fits',
             'part3-2000screen-I=8-2200nm-diffractionLimited.fits',
            'part3-2000screen-I=8-2200nm-update.fits',
            # 'part3-2000screen-I=11_5-2200nm.fits',
            # 'part3-2000screen-I=11_5-2200nm-update.fits'
            ]
            ]

plt.close('all')
psfProftextPos30mas = [[30,1*10**-4],[30,2*10**-2],[30,6*10**-2]]
psfProftextPos90mas = [[90,1*10**-3],[90,3*10**-3],[90,1*10**-2]]
binnedaveragePos30mas = [[30,5*10**-2],[30,5*10**-2],[30,5*10**-2]]
binnedaveragePos90mas = [[90,2*10**-3],[90,3*10**-3],[90,8*10**-5]]
lim1 = 1*10**-2
lim2 = 1.6*10**-3
baseline = [[25,85],[35,95]],[[lim1,lim2],[lim1,lim2]]
lim1 = 1*10**-3
lim2 = 2.0*10**-4
goal = [[25,85],[35,95]],[[lim1,lim2],[lim1,lim2]]

magnitude = ['diffraction limit',8,15,14]
dict2share = {}
# toto = 
for mag in range(1):
    if dict2share != {} and magnitude[mag] in dict2share.keys() :
        dict2share[magnitude[mag]][wl[waveNb]]=[]
    else:
        dict2share[magnitude[mag]] = {wl[waveNb]:[]}

    # with fits.open('part3-2000screen-1000nm.fits') as hdul1:
    with fits.open(fileList[waveNb][mag]) as hdul1:
        # hdul1.info()
        test = hdul1[0].data
    
    test *= 1/np.max(test)
    # plt.figure(0)
    # plt.clf()
    # plt.imshow(np.log(test))
    imgSize = test.shape[0]
    radialStop = 0.1 #arcsec
    binSize = 0.010#arcsec
    nbBin = np.floor(radialStop/binSize).astype(int)
    binEdges = np.arange(binSize/2,radialStop,binSize)
    
    x = np.arange(-imgSize/2,imgSize/2)+0.5
    # pix2arcsec = 1000*10**-9/39 *180*3600/(2*np.pi)
    # pix2arcsec = 0.00125
    pix2arcsec = hdul1[0].header['PIXELSCL']
    x*=pix2arcsec
    nbPixDiam = np.sum(np.abs(x)<radialStop)
    
    X,Y = np.meshgrid(x,x)
    
    R = np.sqrt(X**2+Y**2)
    cR = R[int(imgSize/2-nbPixDiam/2):int(imgSize/2+nbPixDiam/2),
           int(imgSize/2-nbPixDiam/2):int(imgSize/2+nbPixDiam/2)]
    phi = np.arctan2(Y, X)+np.pi
    cphi = phi[int(imgSize/2-nbPixDiam/2):int(imgSize/2+nbPixDiam/2),
           int(imgSize/2-nbPixDiam/2):int(imgSize/2+nbPixDiam/2)]
    fieldStop = cR<radialStop
    cR *= fieldStop
    ctest = test[int(imgSize/2-nbPixDiam/2):int(imgSize/2+nbPixDiam/2),
            int(imgSize/2-nbPixDiam/2):int(imgSize/2+nbPixDiam/2)]
    
    xvar = np.arange(0,binEdges[-1],binSize)*10**3
    
    
    xx,yy = np.meshgrid(np.arange(nbPixDiam),np.arange(nbPixDiam))
    
    popt, pcov = opt.curve_fit(twoD_Gaussian, (xx.ravel(), yy.ravel()), 
                               ctest.ravel(), p0=(np.max(ctest),nbPixDiam/2,
                                               nbPixDiam/2,1,1,0,0))
    #calculate the reference
    mask = R<=binSize/2
    ref = np.sum(test*mask)
    radDist = np.round(radialStop/binSize).astype(int)
    # trace = np.zeros(test.shape)
    mmin = np.zeros(radDist)
    mmax= np.zeros(radDist)
    maverage = np.zeros(radDist)
    mstd = np.zeros(radDist)
    mmin[0] = mmax[0] = maverage[0] = ref
    
    for rd in range(1,radDist):
        alpha = np.arctan(1/(4*rd))
        nbPos = np.round(2*np.pi/alpha).astype(int)
        curring = []
        for ap in range(nbPos):
            
            tmpX,tmpY = np.meshgrid(x-(rd*binSize*np.cos(alpha*ap)),
                                    x-(rd*binSize*np.sin(alpha*ap)))
            tmpR = np.sqrt(tmpX**2+tmpY**2)
            mask=tmpR<=binSize/2
            curring.append(np.sum(test*mask))
        plt.figure(10+rd)
        plt.clf()
        plt.imshow(np.log(test),cmap = 'hot',extent = [-200*pix2arcsec,200*pix2arcsec
                                               ,-200*pix2arcsec,200*pix2arcsec])
        plt.scatter(rd*binSize*np.cos(alpha*np.arange(nbPos))
                    ,rd*binSize*np.sin(alpha*np.arange(nbPos)),
                    c=np.log(curring),cmap='hot')
        plt.title('radial distance {:.0f} mas'.format(radDist*rd))
            # plt.clf()
            # plt.imshow(trace)
            # plt.pause(0.1)
        # plt.figure(10+rd)
        # plt.clf()
        # plt.hist(curring/ref)
        # plt.title('radial distance {:.0f} mas'.format(radDist*rd))
        # plt.xlabel('contrast values')
        # plt.ylabel('number of samples')
        # mmin[rd] = np.min(curring)
        # mmax[rd] = np.max(curring)
        # maverage[rd] = np.mean(curring)
        # mstd[rd] = np.std(curring)
        
        
    
    plt.figure(5+mag)
    plt.clf()
    plt.semilogy(xvar,mmin/ref,drawstyle = 'steps-mid')
    plt.plot(xvar,mmax/ref,drawstyle = 'steps-mid')
    toto = plt.errorbar(xvar, maverage/ref, yerr = mstd/ref,capsize=5,drawstyle = 'steps-mid')
    plt.plot(baseline[0],baseline[1],'r',linewidth = 3)
    plt.plot(goal[0],goal[1],'m',linewidth = 3)
    plt.xlabel('distance from PSF center [mas]')
    plt.ylabel('relative flux')
    
    plt.title('wl = {:.0f} nm, I={}'.format(hdul1[0].header['WAVE0']*10**9,magnitude[mag]))
    custom_lines = [Line2D([0],[0],color='C0'),
                    Line2D([0],[0],color='C1'),
                    Line2D([0],[0],color='C2'),
                    Line2D([0],[0],color = 'r'),
                    Line2D([0],[0],color = 'm')
                    ]
    plt.legend(custom_lines,['min','max','average+std',
                                'TRS Baseline',
                                'GOAL',
                             ],
               handletextpad=0.2,
               labelspacing = 0.2)
    plt.ylim([9*10**-6,1.5])
    plt.grid()
    plt.tight_layout()
    
    # plt.figure(3)
    # plt.clf()
    # plt.imshow(np.log(ctest))
    # fitted = twoD_Gaussian((xx.ravel(),yy.ravel()),*popt)
    # plt.contour(fitted.reshape(nbPixDiam,nbPixDiam))
    
    profile,radDist = computeRadialProfile(ctest, popt[1], popt[2])
    plt.figure(4)
    # plt.clf()
    nProfile = profile/np.max(profile)
    
    
    mas = radDist*pix2arcsec*10**3
    dict2share[magnitude[mag]][wl[waveNb]]=[list(mas),list(nProfile)]
    plt.semilogy(mas,nProfile)
    plt.xlim([0,100])
    # plt.set_yscale('log')
    plt.xlabel('distance from PSF center [mas]')
    plt.ylabel('relative intensity')
    
    # argclosest = np.argmin(np.abs(mas-90))
    # plt.annotate('I={:.0f} ~ {:.01e}'.format(magnitude[mag],nProfile[argclosest]), 
    #              psfProftextPos90mas[mag])
    # argclosest = np.argmin(np.abs(mas-30))
    # plt.annotate('I={:.0f} ~ {:.01e}'.format(magnitude[mag],nProfile[argclosest]),
    #              psfProftextPos30mas[mag])
    plt.grid()
    plt.title('wavelength {:.0f} nm'.format(hdul1[0].header['WAVE0']*10**9))
    plt.tight_layout()
plt.plot(baseline[0],baseline[1],'r',linewidth = 3)
plt.plot(goal[0],goal[1],'m',linewidth = 3)
plt.grid()
custom_lines = [Line2D([0],[0],color='C0'),
                Line2D([0],[0],color='C1'),
                Line2D([0],[0],color = 'r'),
                Line2D([0],[0],color = 'm')
                ]
plt.ylim([10**-5,1.5*10**0])
plt.legend(custom_lines,['diffraction limit','I=8'
                           ,'Baseline','goal'
                         ],handletextpad=0.2,
               labelspacing = 0.2)
# plt.semilogy([90,90], [5*10**-5,1],'k--')
# plt.semilogy([30,30], [5*10**-5,1],'k--')