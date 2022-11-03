# -*- coding: utf-8 -*-
"""
Created on Wed Jul  6 15:45:10 2022

@author: anne.cheffot
"""
import numpy as np
import matplotlib.pyplot as plt
import json
from radialProfile import computeRadialProfile, binnedContrast, plotBinnedContrast
import itertools as itt
from astropy.io import fits

plt.close('all')

wavelength = [1000e-9]#,1600*10**-9,2200*10**-9]
# atmserie = ['I=1','I=2','I=3','I=4','I=5']#,'I=6','I=7']
atmserie = ['I=2','I=5']
# atmserie = ['seeing1']#,'seeing2','seeing3','seeing4','seeing5','seeing6','seeing7','seeing8']
rmsStaticPist = [200e-9]
# rmsdynPist = [25,50,75,100,150,200,300,400]
# rmsdynPist = [np.round(a*1e-9,9) for a in rmsdynPist]
rmsdynPist = [50e-9]
# rmsTipTilt = [25,50,75,100,150,200,300,400]
# rmsTipTilt = [np.round(a*1e-9,9) for a in rmsTipTilt]
rmsTipTilt = [100e-9]
displayAllSeeds = False
savefigFlag = False

pathToData = 'dataPackage/'

maxRadDist = 0.25 #mas

allImg = {}
keyCascade = itt.product(wavelength,atmserie,rmsStaticPist,rmsdynPist,rmsTipTilt)
#prepare the slots for the images
for a in keyCascade:
    allImg.update({a:{'data':[],'pixelScale':[]
                      ,'bmax':[[],[]],'bmin':[[],[]]
                      ,'serialNumber':[],'diffLim':[]}})

fluxFrame = [49826,1251.57,198.36,31.44,4.98,0.79,0.13]
# fluxFrame = [a*2 for a in fluxFrame]
seeingConv = [0.2,0.44,0.65,0.85,1.06,1.2,1.5,2.0]
# fluxFrame = [78.97]*4
zeromoff = 15.74
magnitude = np.round(-2.5*np.log10(fluxFrame)+zeromoff)
Nstar = 10**12
netoile = 0.1
snr = 2

with open(pathToData+'dataPackage.json','r') as file:
    listFits = json.load(file)

for lf in listFits:
    hdul = fits.open(pathToData+lf)
    
    #search the atmosphere serie name used
    try:
        indAtmSerie = [a for a in range(len(atmserie)) 
                       if hdul[0].header['ATMSERIE'].count(atmserie[a])>0][0]
    except KeyError:
        indAtmSerie = [a for a in range(len(atmserie)) 
                       if hdul[0].header['MAGNITUD'].count(atmserie[a])>0][0]
    except:
        continue
    
    try:
        curComb = (np.round(hdul[0].header['WAVELEN'],9)
                   ,atmserie[indAtmSerie]
                   ,np.round(hdul[0].header['RMSSTAT'],9)
                   ,np.round(hdul[0].header['RMSSEGJI'],9)
                   ,np.round(hdul[0].header['RMSWS'],9))
    except KeyError:
        curComb = (np.round(hdul[0].header['WAVELEN'],9)
                   ,atmserie[indAtmSerie]
                   ,np.round(hdul[0].header['RMSSTAT'],9)
                   ,np.round(hdul[0].header['RMSJIT'],9)
                   ,np.round(rmsTipTilt[0],9))
    
    if curComb in allImg.keys():
        allImg[curComb]['data'].append(hdul[0].data)
        allImg[curComb]['pixelScale'].append(hdul[0].header['PIXELSCL'])
        allImg[curComb]['serialNumber'].append(lf.split('_')[0])
        allImg[curComb]['diffLim'].append(hdul[0].header['DIFFLMT'])
    
    del hdul[0].data
    hdul.close()

plt.figure(0)
plt.clf()
plt.figure(1)
plt.clf()

colour = 0

for k in allImg.keys():
    allImg[k]['data'] = np.array(allImg[k]['data'])
    
    if displayAllSeeds:
        plt.figure(0)
        for i in range(allImg[k]['data'].shape[0]):
            prof,xpos = computeRadialProfile(allImg[k]['data'][i]
                                             , allImg[k]['data'][i].shape[0]/2
                                             , allImg[k]['data'][i].shape[1]/2)
            xpos*=allImg[k]['pixelScale'][i]
            prof/=np.max(prof)
            plt.loglog(xpos[xpos<=maxRadDist]*10**3,prof[xpos<=maxRadDist],'C{:.0f}'.format(colour))
            
    # plt.figure(1)
    
    # for i in range(allImg[k]['data'].shape[0]):
    # # for i in range(2):
    #     print(k,i)
    #     cropSize = int(100*10**-3/allImg[k]['pixelScale'][i])
    #     img = allImg[k]['data'][i].copy()
    #     tmpimg = img[int(img.shape[0]/2-cropSize/2):int(img.shape[0]/2+cropSize/2)
    #                   ,int(img.shape[0]/2-cropSize/2):int(img.shape[0]/2+cropSize/2)]
    #     bdata = binnedContrast(tmpimg, 0.01, allImg[k]['pixelScale'][i])
        
    #     allImg[k]['bmin'][0].append(bdata[1])
    #     allImg[k]['bmax'][0].append(bdata[2])
    #     plt.plot(bdata[0],bdata[1],drawstyle = 'steps-mid')
    #     plt.plot(bdata[0],bdata[2],drawstyle = 'steps-mid')
        
    #     cropSize = int(1000*10**-3/allImg[k]['pixelScale'][i])
    #     img = allImg[k]['data'][i].copy()
    #     tmpimg = img[int(img.shape[0]/2-cropSize/2):int(img.shape[0]/2+cropSize/2)
    #                   ,int(img.shape[0]/2-cropSize/2):int(img.shape[0]/2+cropSize/2)]
    #     bdata = binnedContrast(tmpimg, 0.1, allImg[k]['pixelScale'][i])
        
    #     allImg[k]['bmin'][1].append(bdata[1])
    #     allImg[k]['bmax'][1].append(bdata[2])
    
    plt.figure(0)
    if allImg[k]['data'].shape[0]>1:
        avrImg = np.sum(allImg[k]['data'],axis=0)
    else:
        avrImg = allImg[k]['data'][0].copy()
    prof,xpos = computeRadialProfile(avrImg, avrImg.shape[0]/2, avrImg.shape[1]/2)
    xpos*= np.mean(allImg[k]['pixelScale'])
    prof/=np.max(prof)
    # for ld in range(2,5):
    for ld in [2,3,4,10]:
        diffLim = np.mean(allImg[k]['diffLim'])
        tmp = np.argmin(np.abs(xpos/diffLim-ld))
        print(k)
        print(r"radial ditance {:.02f} ({:.02f} $\lambda/D$), contrast: {:.06f}"
              .format(xpos[tmp]*10**3,xpos[tmp]/diffLim,prof[tmp]))
    
    plt.loglog(xpos[xpos<=maxRadDist]*10**3,prof[xpos<=maxRadDist],'C{:.0f}'.format(colour)
                ,label = 'I={}'.format(magnitude[int(k[1][-1])-1])
                # ,label = '{}"'.format(seeingConv[int(k[1][-1])-1])
                # ,label = '{:.0f} nm'.format(k[4]*10**9)
               # ,label='{:.0f} nm'.
               )
    
    plt.figure(1)
    plt.loglog(xpos[xpos<=maxRadDist]*10**3, snr*np.sqrt(prof[xpos<=maxRadDist]/(netoile*Nstar)))
    
    
    cropSize = int(100*10**-3/np.mean(allImg[k]['pixelScale']))
    img = avrImg.copy()
    tmpimg = img[int(img.shape[0]/2-cropSize/2):int(img.shape[0]/2+cropSize/2)
                  ,int(img.shape[0]/2-cropSize/2):int(img.shape[0]/2+cropSize/2)]
    bdata = binnedContrast(tmpimg, 0.01, np.mean(allImg[k]['pixelScale']))
    # bdata = list(bdata)
    # bdata[1]=np.min(allImg[k]['bmin'][0],axis=0)
    # bdata[2] = np.max(allImg[k]['bmax'][0],axis=0)
    
    plotBinnedContrast(bdata, k[0]*10**9, 
                            'I={}'.format(magnitude[int(k[1][-1])-1]),
                            # '{}"'.format(seeingConv[int(k[1][-1])-1]),
                            # 'psf tip tilt {:.0f} nm'.format(k[4]*10**9),
                            figName='10mas_{}_wl{:.0f}nm_{}_st{:.0f}nm_dp{:.0f}nm_tt{:.0f}nm.png'
                            .format(allImg[k]['serialNumber'][0],k[0]*10**9,k[1]
                                    ,k[2]*10**9,k[3]*10**9,k[4]*10**9)
                            ,saveFig=savefigFlag
                            , ylimp=[10**-4,1.1]
                            )
    
    cropSize = int(1000*10**-3/np.mean(allImg[k]['pixelScale']))
    img = avrImg.copy()
    tmpimg = img[int(img.shape[0]/2-cropSize/2):int(img.shape[0]/2+cropSize/2)
                  ,int(img.shape[0]/2-cropSize/2):int(img.shape[0]/2+cropSize/2)]
    bdata = binnedContrast(tmpimg, 0.1, np.mean(allImg[k]['pixelScale']))
    # bdata = list(bdata)
    # bdata[1]=np.min(allImg[k]['bmin'][0],axis=0)
    # bdata[2] = np.max(allImg[k]['bmax'][0],axis=0)
    
    plotBinnedContrast(bdata, k[0]*10**9, 
                            # 'I={}'.format(magnitude[int(k[1][-1])-1]),
                            '{}"'.format(seeingConv[int(k[1][-1])-1]),
                            # 'psf tip tilt {:.0f} nm'.format(k[4]*10**9),
                            figName='100mas_{}_wl{:.0f}nm_{}_st{:.0f}nm_dp{:.0f}nm_tt{:.0f}nm.png'
                            .format(allImg[k]['serialNumber'][0],k[0]*10**9,k[1]
                                    ,k[2]*10**9,k[3]*10**9,k[4]*10**9)
                            ,saveFig=savefigFlag
                            , ylimp=[10**-4,1.1]
                            )
    
    colour +=1




#open the DL data
with open('DLpsf.json','r') as file:
    listDLFile = json.load(file)

dlwl = [int(a.split('_')[-1].split('n')[0]) for a in listDLFile]
plt.figure(0)
for fi in range(len(listDLFile)):
    if dlwl[fi] in [int(a*10**9) for a in  wavelength]:
        hdul = fits.open(listDLFile[fi])
        img = hdul[0].data.copy()
        pixScale = hdul[0].header['PIXELSCL']
    
        #psf profile
        radProf,xpos = computeRadialProfile(img, (img.shape[0]/2-1), (img.shape[1]/2*1))
        xpos*=pixScale
        radProf /= np.max(radProf)
        del hdul[0].data
        hdul.close()
        
        plt.loglog(xpos[xpos<=maxRadDist]*10**3,radProf[xpos<=maxRadDist],'C{:.0f}'.format(colour),label = 'DL')
        
        colour+=1
        #strehl ratio computation
        for k in allImg.keys():
            if allImg[k]['data'].shape[0]>1:
                avrImg = np.sum(allImg[k]['data'],axis=0)/allImg[k]['data'].shape[0]
            else:
                avrImg = allImg[k]['data'][0].copy()
            print('case:{}, SR : {}'.format(k,np.max(avrImg)/np.max(img)))

plt.xlabel('radial distance (mas)')
plt.ylabel('Contrast')
plt.title('wavelength {:.0f} nm'.format(wavelength[0]*10**9))
plt.ylim([10**-8,1.2])
plt.grid()
diffLim = hdul[0].header['DIFFLMT']*10**3
# plt.semilogx([diffLim,diffLim],[10**-8,1],'--')
plt.semilogx([2*diffLim,2*diffLim],[10**-8,1],'k--')
plt.semilogx([4*diffLim,4*diffLim],[10**-8,1],'k--')
plt.semilogx([10*diffLim,10*diffLim],[10**-8,1],'k--')
plt.legend(ncol=2,columnspacing=0.5)
plt.tight_layout()




#%% test sigma 

# X,Y = np.meshgrid(np.linspace(-avrImg.shape[0]/2,avrImg.shape[0]/2,avrImg.shape[0]),
#                   np.linspace(-avrImg.shape[1]/2,avrImg.shape[1]/2,avrImg.shape[1]))
# X*=allImg[k]['pixelScale'][0]/(wavelength[0]/38.54*180*3600/np.pi)
# Y*=allImg[k]['pixelScale'][0]/(wavelength[0]/38.54*180*3600/np.pi)

# R = np.sqrt(X**2+Y**2)
# # plt.close('all') 
# # plt.imshow(R<1)
# test = avrImg/np.sum(avrImg)

# simpleStd = []

# for a in range(1,11):
#     simpleStd.append(np.std(test[np.logical_and(R>=a-1,R<a)]))
    
# plt.figure()
# plt.semilogy(np.arange(1,11),simpleStd)
# listPixel = np.argwhere(np.logical_and(R<11.5,R>=0.5))
# stdHC = []
# for p in listPixel:
#     print('\r {}'.format(p),end = '')
#     Rtmp = np.sqrt((X-X[p[0],p[1]])**2+(Y-Y[p[0],p[1]])**2)
#     incluPix = np.logical_and(R>(R[p[0],p[1]]-0.5),R<=(R[p[0],p[1]]+0.5))
#     excluPix = Rtmp>0.6
#     stdHC.append(np.std(test[np.logical_and(incluPix,excluPix)]))
    

# listPixel = np.array(listPixel)
# plt.semilogy(R[listPixel[:,0],listPixel[:,1]],stdHC,'+')

