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
# magnitude = [a for a in range(1,8)]
# magnitude = [2,5]
magnitude = [1]
# staticPist = [100*10**-9,200*10**-9]
staticPist = [25,50,75,100,150,200,300,400]
# staticPist = [200]
staticPist = [a*10**-9 for a in staticPist]
dynamicPist = 0
psfJitter = 0
# delay = [a for a in range(1,9)]
# wavelength = [1.0*10**-6,1.6*10**-6,2.2*10**-6]
wavelength = [1.0*10**-6]
seed = [12345,47912,35078,30645,4005,34038,91005,16309,19815,7845,89338,17402,
        75086,78130,99707,8679,94826,55305,88253,65789]
nPix = 1600
nScreen = 100

listImag=[]
combination = itt.product(wavelength,magnitude,staticPist,seed)
for comb in combination:
    wl,mag,sp,sd = comb
    print('wavelength : {}nm, magnitude {}'.format(wl*10**9,mag))
    # outFileSuf = 'DP-I={}'.format(mag)
    # outFileSuf = 'DL'
    # outFileSuf = 'statPist{:.0f}'.format(sp*10**9)
    outFileSuf = 'test'
    atmPref = 'resOPDVarFlux/I={}/ResidualOPD_mag_{}_iteration_'.format(mag,mag)
    # atmPref = 'LATENCY/delay{}/ResidualOPD_delay_{}_iteration_'.format(mag,mag)
    # atmPref = 'SEEING/seeing{}/ResidualOPD_seeing_{}_iteration_'.format(mag,mag)
    
    outFile = SCAOSim(wl, nPix, nScreen, sp, dynamicPist, psfJitter, 0.0
                      , sd, ppMap,atmFilePrefix = atmPref
                      ,outputFileSufix=outFileSuf,esoData=False,
                      path2esoData='C:/Users/anne.cheffot/simulations/python/ANDES-SCAO/',
                      esoPsfTipTilt='Time_hist_TT_wind10ms_TelZ45M2_290sec_500Hz.fits')
    
    listImag.append(outFile)

with open('testPsf.json', 'w') as file:
    json.dump(listImag,file)
