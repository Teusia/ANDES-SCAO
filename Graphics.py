# -*- coding: utf-8 -*-
"""
Created on Fri May 27 10:54:46 2022

@author: esoria

Graphics generator
"""
from demo_Scao import *
from ev_PSF import *
import matplotlib.pyplot as plt
import pandas as pd

#Paramaters:
I=8
Band= 'Y band'


#Requirements for each WL
data = {
  "Wavelength": [1*10**-6, 1.6*10**-6, 2.2*10**-6],
  "Baseline top": [0.01, 0.01, 0],
  "Baseline bottom": [0.0017, 0.001, 0],
  "Goal top": [0.001, 0.001,0],
  "Goal bottom": [0.0002, 0.0001, 0],
}


df = pd.DataFrame(data, index = ["Y band", "H band", "K band"])

 
wavelength= df.at[Band,'Wavelength']
base_top=df.at[Band,'Baseline top']
base_bot=df.at[Band,'Baseline bottom']
goal_top=df.at[Band,'Goal top']
goal_bot=df.at[Band,'Goal bottom']




#COMPARE THE DIFFERENT CONTRIBUTIONS
DL,rms_DL=SCAOSim(wavelength = wavelength , nPix = 400, nScreen = 2, SDIslandEffect = 0, SDSegmentJitter = 0
            ,rmsWindshake=0
            , mag=0, pupilFileName='Tel-Pupil.fits', atmFilePrefix='screen')

pixDL,radDL=ev_PSF(DL)

WS,rms_WS=SCAOSim(wavelength = wavelength , nPix = 400, nScreen = 2, SDIslandEffect = 0, SDSegmentJitter = 0
           ,rmsWindshake=1e-07
           ,mag = I , pupilFileName='Tel-Pupil.fits', atmFilePrefix='screen')
pixWS,radWS=ev_PSF(WS)

atm,rms_atm=SCAOSim(wavelength = wavelength , nPix = 400, nScreen = 2, SDIslandEffect = 0, SDSegmentJitter = 0
            ,rmsWindshake=0
            ,mag = I, pupilFileName='Tel-Pupil.fits', atmFilePrefix='screen')
pix,rad=ev_PSF(atm)
atm_p,rms_atm_p=SCAOSim(wavelength = wavelength , nPix = 400, nScreen = 2, SDIslandEffect = 0.15, SDSegmentJitter = 0
            ,rmsWindshake = 0
            ,mag = I, pupilFileName='Tel-Pupil.fits', atmFilePrefix='screen')
pix1,rad1=ev_PSF(atm_p)


atm_v,rms_atm_v=SCAOSim(wavelength = wavelength , nPix = 400, nScreen = 2, SDIslandEffect = 0, SDSegmentJitter = 0.05
            ,rmsWindshake = 0
            ,mag = I, pupilFileName='Tel-Pupil.fits', atmFilePrefix='screen')
pix2,rad2=ev_PSF(atm_v)

atm_p_v,rms_atm_p_v=SCAOSim(wavelength = wavelength , nPix = 400, nScreen = 2, SDIslandEffect = 0.15, SDSegmentJitter = 0.05
            ,rmsWindshake = 0
            ,mag = I, pupilFileName='Tel-Pupil.fits', atmFilePrefix='screen')
pix3,rad3=ev_PSF(atm_p_v)
atm_p_v_w,rms_atm_p_v_w=SCAOSim(wavelength = wavelength , nPix = 400, nScreen = 2, SDIslandEffect = 0.15, SDSegmentJitter = 0.05
            ,rmsWindshake = 1e-07
            ,mag = I, pupilFileName='Tel-Pupil.fits', atmFilePrefix='screen')
pix4,rad4=ev_PSF(atm_p_v_w)

#plot
plt.figure()
plt.yscale("log")
#plt.title('Radial distribution, I={}, WL= {}um'.format(I, round(wavelength*10**6,2)), fontsize=15)
plt.xlabel('Distance from PSF center (mas)', fontsize=12)
plt.ylabel('Frequency', fontsize=12)
plt.plot(pixDL[0:80], radDL[0:80],  label='DL')
plt.plot(pixWS[0:80], radWS[0:80],  label='AOres + Ws')

plt.plot(pix[0:80], rad[0:80], label='AOres')
plt.plot(pix1[0:80], rad1[0:80], label='AOres+ static piston')
plt.plot(pix2[0:80], rad2[0:80], label='AOres + dynamic piston')
plt.plot(pix3[0:80], rad3[0:80], label='AOres + static + dynamic piston')
plt.plot(pix4[0:80], rad4[0:80], label='AOres + WS + static + dynamic piston')

plt.hlines(y=base_top, xmin=25 ,xmax=35, colors='aqua', linestyles='--', lw=2, label='Baseline')
plt.hlines(y=base_bot, xmin=85, xmax=95, colors='aqua', linestyles='--', lw=2)
plt.hlines(y=goal_top, xmin=25 ,xmax=35, colors='darkviolet', linestyles='--', lw=2, label='Goal')
plt.hlines(y=goal_bot, xmin=85, xmax=95, colors='darkviolet', linestyles='--', lw=2)
plt.legend(fontsize=8)
plt.grid()
plt.show()



#COMPARE THE INFLUENCE OF THE REFECENCE STAR
plt.figure()
step=0
for i in [0,8,12,14]:
    WS,rms_WS=SCAOSim(wavelength = wavelength , nPix = 400, nScreen = 2, rmsIslandEffect = 0, rmsSegmentJitter = 0
            ,rmsWindshake=0
            ,mag = i , pupilFileName='Tel-Pupil.fits', atmFilePrefix='screen'
            ,outputFileSufix='I={}-{}nm.fits'.format(i,wavelength))
    pixWS,radWS=ev_PSF(WS)
    
    plt.yscale("log")
    plt.grid(True, which="major")
    plt.title('Radial distribution, WL= {}um'.format(round(wavelength*10**6,2)), fontsize=15)
    plt.xlabel('Distance from PSF center (mas)', fontsize=12)
    plt.ylabel('Frequency', fontsize=12)   
    if i ==0:
       
        plt.legend(fontsize=8)
        plt.plot(pixWS[0:140], radWS[0:140],  label='DL')
        plt.hlines(y=base_top, xmin=25 ,xmax=35, colors='aqua', linestyles='--', lw=2, label='Baseline')
        plt.hlines(y=base_bot, xmin=85, xmax=95, colors='aqua', linestyles='--', lw=2)
        plt.hlines(y=goal_top, xmin=25 ,xmax=35, colors='darkviolet', linestyles='--', lw=2, label='Goal')
        plt.hlines(y=goal_bot, xmin=85, xmax=95, colors='darkviolet', linestyles='--', lw=2)
    else:
        plt.plot(pixWS[0:150], radWS[0:150],  label='Mag: {}'.format(i))
        plt.axvline(30, color='black', ls="dotted")
        plt.axvline(90, color='black', ls="dotted")
        step+=20
    plt.annotate('I={}: {}'.format(i,format(radWS[30],'.1E')),
     xy=(0, 0),
     xycoords='figure points',
     xytext=(128,175+step),
     fontsize=11)
    plt.annotate('I={}: {}'.format(i,format(radWS[90],'.1E')),
     xy=(0, 0),
     xycoords='figure points',
     xytext=(225,175+step),
     fontsize=11)
    plt.legend(fontsize=8)
plt.show()


#COMPARE THE INFLUENCE OF THE DISTRUBUTION OF THE SEGMENTS SAME SD
import demo_Scao_tmp
plt.figure()
for i in range(10):
    plt.yscale("log") 
    if i ==0:    
        WS,SD=demo_Scao_tmp.SCAOSim(wavelength = wavelength , nPix = 400, nScreen = 100,rmsIslandEffect=0, rmsSegmentJitter = 0
            ,rmsWindshake=0
            ,mag = 0 , pupilFileName='Tel-Pupil.fits', atmFilePrefix='screen')
        pixWS,radWS=ev_PSF(WS)
   
        plt.plot(pixWS[0:140], radWS[0:140],  label='DL')
        plt.hlines(y=base_top, xmin=25 ,xmax=35, colors='aqua', linestyles='--', lw=2, label='Baseline')
        plt.hlines(y=base_bot, xmin=85, xmax=95, colors='aqua', linestyles='--', lw=2)
        plt.hlines(y=goal_top, xmin=25 ,xmax=35, colors='darkviolet', linestyles='--', lw=2, label='Goal')
        plt.hlines(y=goal_bot, xmin=85, xmax=95, colors='darkviolet', linestyles='--', lw=2)
    else:
        WS,SD=demo_Scao_tmp.SCAOSim(wavelength = wavelength , nPix = 400, nScreen = 100,rmsIslandEffect=0.15, rmsSegmentJitter = 0
            ,rmsWindshake=0
            ,mag = I , pupilFileName='Tel-Pupil.fits', atmFilePrefix='screen')
        pixWS,radWS=ev_PSF(WS)
   
        plt.yscale("log") 
        plt.plot(pixWS[0:140], radWS[0:140])
    
    plt.grid(True, which="major")
    #plt.title('Radial distribution,I:{}  WL= {}um, SD:{} nm'.format(I,round(wavelength*10**6,2),round(SD*1000)), fontsize=15)
    plt.xlabel('Distance from PSF center (mas)', fontsize=12)
    plt.ylabel('Frequency', fontsize=12)
    plt.legend()
    print(i)
plt.show()


#****USEFULL******

#DELATE THE PUPIL FITS FILES
import shutil
import os
dirpath=('C:/Users/esoria.DOMAINT/Desktop/ANDES/vibra/')
shutil.rmtree(dirpath)
os.mkdir(dirpath)

dirpath=('C:/Users/esoria.DOMAINT/Desktop/ANDES/static/')
shutil.rmtree(dirpath)
os.mkdir(dirpath)