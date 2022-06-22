# -*- coding: utf-8 -*-

"""
Created on Wed Jun  1 13:12:05 2022

@author: esoria
"""

from multiprocessing import Pool
from scao_multi_v1 import *
from ev_PSF import *
import matplotlib.pyplot as plt
import numpy as np
global str
import pandas as pd
from functools import partial

#PARAMETERS

target=[0, 1, 1.5, 2, 2.5,3, 3.5] #VALUES SD PISTON STATIC

I= 8 #REF MAG

Band= 'K band'


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

size=np.size(target)
log=np.zeros((size,2))
if __name__ == '__main__':
    p = Pool(4)
    log=p.map(partial(SCAOSim_mul, I, wavelength),target)
    #STARMAP

for i in range(size):    
    j= "u" + str(i) 
    locals()[j] =log[i]
 

plt.figure()
if u0[0] != 0 and u4[1]!=0:
    for i in range(size):    
        j,k = ev_PSF(locals()["u" + str(i)][0]) 
        plt.yscale("log")
        plt.title('Radial distribution, I={}, WL= {}um'.format(I, round(wavelength*10**6,2)), fontsize=15)
        if locals()["u" + str(i)][1]==0:
            plt.plot(j[0:80], k[0:80], label='psf Res AO')
            plt.hlines(y=base_top, xmin=25 ,xmax=35, colors='aqua', linestyles='--', lw=2, label='Baseline')
            plt.hlines(y=base_bot, xmin=85, xmax=95, colors='aqua', linestyles='--', lw=2)
            plt.hlines(y=goal_top, xmin=25 ,xmax=35, colors='darkviolet', linestyles='--', lw=2, label='Goal')
            plt.hlines(y=goal_bot, xmin=85, xmax=95, colors='darkviolet', linestyles='--', lw=2)
        else:
            plt.plot(j[0:80], k[0:80], label='psf + petaling :{}'.format(round(locals()["u" + str(i)][1]*1000,-1)) +str(' nm'))
        
        plt.xlabel('Distance from PSF center (mas)', fontsize=12)
        plt.ylabel('Frequency', fontsize=12)
        plt.grid(True, which="major")
        plt.legend(fontsize=10)
        
    plt.show()    
  

       