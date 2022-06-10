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
import time
from functools import partial
target=[0, 0.2, 0.3, 0.35, 0.4, 0.5, 0.75, 1, 2]
size=np.size(target)
magn=[8,11.5,14,15]
atmWF=[1,1.6,2.16,2.57]
idx = 0
I=magn[idx]
atm=atmWF[idx]
wavelength= 1*10**-6
log=np.zeros((size,2))


if __name__ == '__main__':
    p = Pool(5)
    log=p.map(partial(SCAOSim, atm, wavelength),target)
    #STARMAP

for i in range(size):    
    j= "u" + str(i) 
    locals()[j] =log[i]

plt.figure()
if u0[0] != 0 and u6[1]!=0:
    for i in range(size):    
        j,k = ev_PSF(locals()["u" + str(i)][0]) 
        plt.yscale("log")
        plt.title('Radial distribution, I={}, WL= {}um'.format(I, round(wavelength*10**6,2)), fontsize=15)
        if locals()["u" + str(i)][1]==0:
            plt.plot(j[0:80], k[0:80], label='psf DL')
        else:
            plt.plot(j[0:80], k[0:80], label='psf + petaling RMS:{}'.format(round(locals()["u" + str(i)][1]*1000,-1)) +str(' nm'))
        plt.xlabel('Distance from PSF center (mas)', fontsize=12)
        plt.ylabel('Frequency', fontsize=12)
        plt.legend(fontsize=7)
        plt.grid()
    plt.show()

