# -*- coding: utf-8 -*-
"""
Created on Wed Jul  6 11:05:52 2022

@author: anne.cheffot
"""

import numpy as np
import poppy as po
import datetime
from pupil_ANDES import buildPupil
from demoScao import SCAOSim
import itertools as itt
import json

ppMap = buildPupil()
magnitude = [a for a in range(1,8)]
# magnitude = [2,5]
# magnitude = [1,2,3,4]
# staticPist = [100*10**-9,200*10**-9]
# staticPist = [25,50,75,100,150,200,300,400]
staticPist = [0e-9]
# staticPist = [a*10**-9 for a in staticPist]
# dynamicPist = [25,50,75,100,150,200,300,400]
# dynamicPist= [a*10**-9 for a in dynamicPist]
dynamicPist = 0e-9
# psfJitter =[25,50,75,100,150,200,300,400]
# psfJitter= [a*10**-9 for a in psfJitter]
psfJitter = 0e-9
# delay = [a for a in range(1,9)]
# wavelength = [1.0e-6,1.6e-6,2.2e-6]
# wavelength = [900e-9,800e-9,700e-9,600e-9]
wavelength = [0.75*10**-6]
seed = [12345]#,47912,35078,30645,4005,34038,91005,16309,19815,7845,89338,17402,
            # 75086,78130,99707,8679,94826,55305,88253,65789]
nPix = 1600
nScreen = 2000

listImag=[]
combination = itt.product(wavelength,magnitude,staticPist,seed)
for comb in combination:
    wl,mag,sp,sd = comb
    print(comb)
    # print('wavelength : {}nm, magnitude {}'.format(wl*10**9,mag))
    # outFileSuf = 'DP-I={}'.format(mag)
    # outFileSuf = 'DL'
    # outFileSuf = 'DLvisible{:03.0f}'.format(sp*10**9)
    outFileSuf = 'EPVisi750varflux'
    atmPref = 'FLUX/I={}/ResidualOPD_mag_{}_iteration_'.format(mag,mag)
    # atmPref = 'LATENCY/delay{}/ResidualOPD_delay_{}_iteration_'.format(mag,mag)
    # atmPref = 'SEEING/seeing{}/ResidualOPD_seeing_{}_iteration_'.format(mag,mag)
    
    outFile = SCAOSim(wl, nPix, nScreen, sp, dynamicPist, psfJitter, 1.0
                      , sd, ppMap,atmFilePrefix = atmPref
                      ,outputFileSufix=outFileSuf,esoData=True,
                      path2esoData='/home/alcheffot/ANDES/ANDES-SCAO/',
                      esoPsfTipTilt='Time_hist_TT_wind10ms_TelZ45M2_290sec_500Hz.fits')
    
    listImag.append(outFile)

with open('EPvisi750varFlux.json', 'w') as file:
    json.dump(listImag,file)
