# -*- coding: utf-8 -*-
"""
Created on Wed Jul  6 15:45:10 2022

@author: anne.cheffot
"""
import numpy as np
import matplotlib.pyplot as plt
import json
from radialProfile import computeRadialProfile, binnedContrast, plotBinnedContrast
from astropy.io import fits

plt.close('all')

#open the DL data
with open('DLpsf.json','r') as file:
    listDLFile = json.load(file)

for fi in listDLFile:
    hdul = fits.open(fi)
    img = hdul[0].data.copy()
    pixScale = hdul[0].header['PIXELSCL']
    
    #psf profile
    radProf,xpos = computeRadialProfile(img, (img.shape[0]/2)-1, (img.shape[1]/2)-1)
    xpos*=pixScale*10**3
    radProf *= 1/np.max(radProf)
    hdul.close()
    
    argkeep = xpos<(img.shape[0]/2*pixScale*10**3)
    xpos = xpos[argkeep]
    radProf = radProf[argkeep]
    if '1000nm.' in fi:
        plt.figure(1000)
        plt.loglog(xpos,radProf,label = 'DL')
        # plt.semilogy(xpos,radProf,label = 'I={:.0f}'.format(magnitude[curMag-1]))
    elif '1600nm.' in fi:
        plt.figure(1600)
        plt.loglog(xpos,radProf,label = 'DL')
    elif '2200nm.' in fi:
        plt.figure(2200)
        plt.loglog(xpos,radProf,label = 'DL')
    
    # plt.figure(fi)
    # plt.clf()
    # plt.imshow(np.log(img/np.max(img)),cmap = 'gist_heat',vmin=-15)
    # plt.axis('off')
    # plt.savefig(fi.split('.')[0]+'.svg')

with open('testPsf.json','r') as file:
    listFits = json.load(file)

fluxFrame = [49826,1251.57,198.36,31.44,4.98,0.79,0.13]
# fluxFrame = [78.97]*4
zeromoff = 15.74
magnitude = np.round(-2.5*np.log10(fluxFrame)+zeromoff)
delay = [a for a in range(1,9)]

avgimg = {1000:{},1600:{}, 2200:{},'pixscl':{}}

for fi in listFits :
    # curMag = int(fi.split('_')[1][-1])
    # curMag = int(fi.split('_')[1][-1])
    # curMag = int(fi.split('_')[1][-3:])
    curMag = fi.split('_')[1]
    hdul = fits.open(fi)
    # if hdul[0].header['RMSSTAT']!=100e-9:
    #     continue
    # print(fi)
    wl = hdul[0].header['WAVELEN']*10**9
    img = hdul[0].data.copy()
    pixScale = hdul[0].header['PIXELSCL']
    radProf,xpos = computeRadialProfile(img, (img.shape[0]/2)-1, (img.shape[1]/2)-1)
    xpos*=pixScale*10**3
    radProf *= 1/np.max(radProf)
    hdul.close()
    
    #crop the statitically unwanted data points
    argkeep = xpos<(img.shape[0]/2*pixScale*10**3)
    xpos = xpos[argkeep]
    radProf = radProf[argkeep]
    
    if '1000nm.' in fi:
        plt.figure(1000)
        # plt.loglog(xpos,radProf,
        #            # label = 'I={:.0f}'.format(magnitude[curMag-1])
        #            # label = '{:.0f} frame'.format(delay[curMag-1])
        #            # label = '{:.0f} "'.format(delay[curMag-1])
        #            # label = '{:.0f} nm'.format(curMag)
        #            )
        
        if curMag in avgimg[1000].keys():
            avgimg[1000][curMag]+=img.copy()
        else:
            avgimg[1000][curMag] = img.copy()
            avgimg['pixscl'].update({1000:pixScale})
        
        # plt.semilogy(xpos,radProf,label = 'I={:.0f}'.format(magnitude[curMag-1]))
    elif '1600nm.' in fi:
        plt.figure(1600)
        # plt.loglog(xpos,radProf,
        #             # label = 'I={:.0f}'.format(magnitude[curMag-1])
        #            # label = '{:.0f} frame'.format(delay[curMag-1])
        #            # label = '{:.0f} "'.format(delay[curMag-1])
        #            # label = '{:.0f} nm'.format(curMag)
        #            )
        if curMag in avgimg[1000].keys():
            avgimg[1000][curMag]+=img.copy()
        else:
            avgimg[1000][curMag] = img.copy()
            avgimg['pixscl'].update({1000:pixScale})
        
    elif '2200nm.' in fi:
        plt.figure(2200)
        # plt.loglog(xpos,radProf,
        #             # label = 'I={:.0f}'.format(magnitude[curMag-1])
        #            # label = '{:.0f} frame'.format(delay[curMag-1])
        #            # label = '{:.0f} "'.format(delay[curMag-1])
        #            # label = '{:.0f} nm'.format(curMag)
        #            )
        
        if curMag in avgimg[1000].keys():
            avgimg[1000][curMag]+=img.copy()
        else:
            avgimg[1000][curMag] = img.copy()
            avgimg['pixscl'].update({1000:pixScale})
    
    # if fi.split('_')[1] in ['DP-I=2','DP-I=5']:
    # # if fi.split('_')[1] in ['delay=2','delay=3','delay=4']:
    #     #cropped figure to only select the usefull field
    #     cropSize = int(100*10**-3/pixScale)
    #     tmpimg = img[int(img.shape[0]/2-cropSize/2):int(img.shape[0]/2+cropSize/2)
    #                   ,int(img.shape[0]/2-cropSize/2):int(img.shape[0]/2+cropSize/2)]
    #     binnedData = binnedContrast(tmpimg, 0.010, pixScale)
    #     plotBinnedContrast(binnedData, wl, 
    #                         'I={:.0f}'.format(magnitude[curMag-1]),
    #                         # 'delay={:.0f}'.format(delay[curMag-1]),
    #                         figName='10mas'+fi.split('.')[0]+'.png',saveFig = True
    #                         , ylimp=[10**-4,1.1]
    #                         )
        
    #     #cropped figure to only select the usefull field
    #     cropSize = int(1000*10**-3/pixScale)
    #     tmpimg = img[int(img.shape[0]/2-cropSize/2):int(img.shape[0]/2+cropSize/2)
    #                   ,int(img.shape[0]/2-cropSize/2):int(img.shape[0]/2+cropSize/2)]
    #     binnedData = binnedContrast(tmpimg, 0.10, pixScale)
    #     plotBinnedContrast(binnedData, wl, 
    #                         'I={:.0f}'.format(magnitude[curMag-1]),
    #                         # 'delay={:.0f}'.format(delay[curMag-1]),
    #                         figName='100mas'+fi.split('.')[0]+'.png',saveFig=True
    #                         , ylimp=[10**-4,1.1]
    #                         )
    
    # plt.figure(fi)
    # plt.clf()
    # plt.imshow(np.log(img/np.max(img)),cmap = 'gist_heat',vmin = -15)
    # plt.axis('off')
    # plt.savefig(fi.split('.')[0]+'.svg')


for fifi in [1000,1600,2200]:
    plt.figure(fifi)
    for k in avgimg[fifi].keys():
        radProf,xpos = computeRadialProfile(avgimg[fifi][k], (img.shape[0]/2)-1, (img.shape[1]/2)-1)
        xpos*=avgimg['pixscl'][fifi]*10**3
        radProf *= 1/np.max(radProf)
        
        #crop the statitically unwanted data points
        argkeep = xpos<(img.shape[0]/2*avgimg['pixscl'][fifi]*10**3)
        xpos = xpos[argkeep]
        radProf = radProf[argkeep]
        plt.plot(xpos,radProf,'b--')
    
    
    plt.title('wavelength {} nm'.format(fifi))
    plt.xlabel('radial distance (mas)')
    plt.ylabel('Contrast')
    plt.ylim([10**-8,1.5])
    plt.grid()
    plt.legend(ncol = 3, columnspacing = 0.5, handletextpad=0.2)
    plt.tight_layout()

