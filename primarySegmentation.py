# -*- coding: utf-8 -*-
"""
Created on Mon Nov 14 09:58:58 2022

@author: anne.cheffot
"""
import numpy as np
import matplotlib.pyplot as plt
import poppy
from astropy.io import fits

def createPrimarySegmentation():
    seglist = [a for a in range((3*4**2)+(3*4)+1,(3*15**2)+(3*15)+1)]
    seglist.extend([a for a in range(723,736)])
    seglist.extend([a for a in range(739,752)])
    seglist.extend([a for a in range(755,768)])
    seglist.extend([a for a in range(771,784)])
    seglist.extend([a for a in range(787,800)])
    seglist.extend([a for a in range(803,816)])
    
    seglist.extend([a for a in range(821,831)])
    seglist.extend([a for a in range(838,848)])
    seglist.extend([a for a in range(855,865)])
    seglist.extend([a for a in range(872,882)])
    seglist.extend([a for a in range(889,899)])
    seglist.extend([a for a in range(906,916)])
    
    ap = poppy.dms.HexSegmentedDeformableMirror(rings = 17, flattoflat = 1.2388,gap=0.000,
                                    segmentlist = seglist,
                                  rotation = 0)
    
    # for s in ap.segmentlist:
    #     ap.set_actuator(s,s*10**-9,0,0)
    
    # # ap.set_actuator(ap.segmentlist[0],-1,0,0)
    
    # plt.figure(1)
    # plt.clf()
    npix = 400
    npixfull = np.ceil(npix*ap.rings/(ap.rings-2)-4).astype(int)
    if npixfull%2!=0:
        npixfull+=1
    
    ap.display(npix =npixfull)
    dpix = ((npixfull-npix)/2).astype(int)
    
    toCrop = ap.transmission.copy()
    ap.transmission=toCrop[dpix:dpix+npix,dpix:dpix+npix]
    ap._transmission = ap.transmission>0
    ap.opd = np.zeros(ap.transmission.shape)
    # plt.figure(2)
    # plt.clf()
    # plt.imshow(ap.transmission)
    # plt.figure(3)
    # plt.clf()
    # plt.imshow(ap.opd)
    
    # plt.figure(4)
    # plt.clf()
    # plt.imshow(ap.transmission-ap.opd*10**9)
    
    wavelength = 850*10**-9
    telDia = ap.pupil_diam.value
    
    hdu = fits.PrimaryHDU(ap.transmission)
    hdu.header['WAVELEN'] = wavelength
    hdu.header['DIAM'] = telDia
    hdu.header['DIFFLMT'] = wavelength/telDia
    hdu.header['OVERSAMP'] = 1
    hdu.header['DET_SAMP'] = 1
    hdu.header['PIXELSCL'] = ap._last_pixelscale.value#telDia/ap.transmission.shape[0] 
    hdu.header['PIXUNIT'] = 'meters'
    hdu.header['BUNIT'] = 'meters'
    hdu.header['side'] = ap.side.value
    hdu.header['flatofla']=ap.flattoflat.value
    hdu.header['rings']=ap.rings
    hdu.header['gap'] = ap.gap.value
    hdu.header['oversamp'] = ap.oversample
    hdu.header['ispadded'] = ap.ispadded
    hdu.header['rotation'] = ap.rotation
    hdu.header['segspaci'] = ap._segment_spacing
    hdu.header['lastnpix'] = npix
    hdu.header['lastpixs'] = ap._last_pixelscale.value
    
    hdu1 = fits.ImageHDU(ap.segmentlist)
    hdu1.header['CONTAINS'] = ' segmentlist'
    
    # hdu2 = fits.ImageHDU(ap._transmission)
    # hdu2.header['CONTAINS'] = '_transmission'
    
    hdul = fits.HDUList([hdu,hdu1])
    
    hdul.writeto('primaryMap.fits',overwrite=True)
    
    return ap

class ELTprimary(poppy.AnalyticOpticalElement):
    def __init__(self,fitsPupilMap, *args, **kwargs):
        """ If your optic has adjustable parameters, then save them as attributes here """
        super().__init__(**kwargs)
        with fits.open(fitsPupilMap,memmap=False) as hdul:
            # hdul = fits.open(fitsPupilMap)
            self.transmission = hdul[0].data.copy()
            del hdul[0].data
            self.pupil_diam = hdul[0].header['DIAM']
            self.npix=hdul[0].header['lastnpix']
            self._last_npix=hdul[0].header['lastnpix']
            self.pixelscale = hdul[0].header['lastpixs']*poppy.u.m/poppy.u.pix
            self._last_pixelscale= hdul[0].header['lastpixs']
            self.rings = hdul[0].header['rings']
            self.side = hdul[0].header['side']*poppy.u.m
            self.gap = hdul[0].header['gap']
            self.flattoflat = hdul[0].header['flatofla']
            
            
            
            self.segmentlist = hdul[1].data
            del hdul[1].data
            del hdul
        
        self._surface = np.zeros((self._n_aper_inside_ring(self.rings + 1)+1, 3))
        self._seg_mask = self.transmission.copy()
        self._transmission = (self.transmission>0).astype(float)
        self._seg_indices = dict()
        for i in self.segmentlist:
            wseg = np.where(self._seg_mask == i+1)
            self._seg_indices[i] = wseg
        
    
    def _n_aper_in_ring(self, n):
        """ How many hexagons or circles in ring N? """
        return 1 if n == 0 else 6 * n

    def _n_aper_inside_ring(self, n):
        """ How many hexagons or circles interior to ring N, not counting N?"""
        return sum([self._n_aper_in_ring(i) for i in range(n)])
    
    
    def set_segment(self,segnum,segmentMove):

        if segnum not in self.segmentlist:
            raise ValueError("Segment {} is not present for this DM instance.".format(segnum))
        self._surface[segnum,0] = segmentMove
        

    def get_opd(self,wave):
        # print(wave.__dict__)
        y, x = self.get_coordinates(wave)
        
        self.opd = np.zeros(self.transmission.shape)
        for s in self.segmentlist:
            self.opd[self.transmission==s]= self._surface[s,0]
        # opd = some_function(x,y, wave.wavelength, self)
        return self.opd

    def get_transmission(self, wave):
        y, x = self.get_coordinates(wave)
        self._transmission = (self.transmission >0).astype(int)
        return self._transmission

    # behind the scenes poppy  will calculate:
    #    phasor = transmission = np.exp(1.j * 2 * np.pi / wave.wavelength * opd)


if __name__=='__main__':
    testap = createPrimarySegmentation()
    prim = ELTprimary('primaryMap.fits')
    plt.clf()
    plt.imshow(testap.transmission)
    plt.show()
    
    plt.figure(2)
    plt.clf()
    rng = np.random.default_rng()
    
    for ii in range(10):
        pist = rng.normal(0,1,len(prim.segmentlist))*100*10**-9
        for i in range(len(prim.segmentlist)):
            prim.set_segment(prim.segmentlist[i],pist[i])
        
        wave = poppy.Wavefront(npix = 400)
        # prim.display(what='opd')
        prim.get_opd(wave)
        plt.imshow(prim.opd)
    
    pupilMask = fits.open('Tel-Pupil300.fits')
    
    plt.figure(3)
    plt.clf()
    plt.imshow(prim._transmission+pupilMask[0].data)
    
    




