# -*- coding: utf-8 -*-
"""
Created on Mon May 30 09:58:19 2022

@author: anne.cheffot
aimed at developping the parallelisation script
"""
import numpy as np
import matplotlib.pyplot as plt
from pupil_ANDES import buildPupil,Mod_segment
from demoScao import SCAOSim
import itertools as itt
import os
import shutil
import datetime

ppMap = buildPupil()

wavelength = [1*10**-6,1.6*10**-6,2.2*10**-6]
seed = [12345]
worseFact = [1.0]
atmPref = ['screen_I=8/screenI=8','screen_I=12/screenI=12']

allconf = itt.product(wavelength,seed,worseFact,atmPref)
worseFact = [1.6]
atmPref = ['screen_I=12/screenI=12']
allconf = itt.chain(allconf,itt.product(wavelength,seed,worseFact,atmPref))

worseFact = [0]
allconf = itt.chain(allconf,itt.product(wavelength,seed,worseFact,atmPref))

nPix = 400
nscreen = 2000
rmsstat = 500*10**-9
rmsjit = 50*10**-9
rmsshake = 0

outList = []

folderName = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
os.mkdir(folderName)

for wl,sd,wf,ap in allconf:
    #create the outFileSuffix depending on the atmScreenPrefix and the worsening factor
    if wf ==0:
        outFileS = 'DL'
    elif 'I=8' in ap:
        outFileS = 'I=8'
    elif ('I=12' in ap) and (wf ==1.0):
        outFileS = 'I=12'
    elif ('I=12' in ap) and (wf ==1.6):
        outFileS = 'I=14'
    
    
    outFile = SCAOSim(wl, nPix, nscreen, rmsstat, rmsjit, rmsshake, wf
                      , sd, ppMap,atmFilePrefix=ap,outputFileSufix=outFileS)
    
    shutil.move(outFile,folderName+'/'+outFile)
    outList.append(outFile)
    

#write in a text file what was done
with open(folderName+'/params.txt','w') as file:
    file.write('list of parameter for the simulation stored in this folder\n')
    file.write('static parameters\n\n')
    file.write( 'nPix = {} \n'.format(nPix))
    file.write( 'nscreen = {}\n'.format(nscreen))
    file.write( 'rmsstat = {}[m]\n'.format(rmsstat))
    file.write( 'rmsjit = {}[m]\n'.format(rmsjit))
    file.write( 'rmsshake = {}[m]\n'.format(rmsshake))
    
    file.write('varied parameters')
    file.write( 'wavelength = {} [m] \n'.format(wavelength))
    file.write( 'seed = {} \n'.format(seed))
    file.write( 'worseFact = {} \n'.format(worseFact))
    file.write( 'atmPref = {}\n'.format(atmPref))
    file.write('\n')
    
    






