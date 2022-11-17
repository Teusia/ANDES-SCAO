"""
Created on Fri Jan 28 11:45:35 2022

@author: anne.cheffot
@author: EstherSoria

"""
import numpy as np
import poppy as po
import os
import datetime
from pupil_ANDES import *
from astropy.io import fits
from primarySegmentation import ELTprimary

def SCAOSim(wavelength,nPix,nScreen,rmsIslandEffect,rmsSegmentJitter,rmsWindshake
            ,atmWorseningFactor,seed,pupilMap,pupilFileName='Tel-Pupil.fits', 
            atmFilePrefix='screen',outputFileSufix='I=8',esoData=False,path2esoData=None,
            esoPsfTipTilt=None,complexOutWF = False,fieldOfView=0.5,pixelScale=0.00125,
            saveIntermediateComplexWavefront=False, pathFileComplex=None,
            primaryFits = None,primaryPistRms = 0):
    '''
    Parameters
    ----------
    wavelength : float
        wavelength of the final image in meters
    nPix : int
        number of pixels across one image
    nScreen : int
        number of time light is propagated through an aberated screen.
    rmsIslandEffect : float
        root mean square error of the island effect in meters
    rmsSegmentJitter : float
        root mean square of the segment jitter in meters
    rmsWindshake : float
        root mean square of the tip tilt induced by the winshake in meters
    atmWorseningFactor : float
        a number to change the atmospheric residuals amplitude (no unit)
    seed : float
        seed for the random noise generator
    pupilFileName : string, optional
        The name of the file in which the pupil shape is stored. The default is 
        'Tel-Pupil.fits'.
    atmFilePrefix : string, optional
        The atmosheric screens should have a prefix, for example 'screen' and be 
        followed by an iterger. The default is 'screen'.
    outputFileSufix : string, optional
        The output file should be date_time followed by a suffix.
        The default is 'I=8-2200nm.fits'.
    esoData : boolean, optionnal
        whether to use the tip tilt informations provided by ESO or not, default : False
    path2esoData : string, optionnal
        Path to where is the eso tiptilt data stored. default: None
    fieldOfView: float,optionnal
        field of view diameter of the image in arc seconds. default: 0.5
        note: fieldOfView/pixelScale must be the number of pixels accross the pupil file(400)
    pixelScale: float, optionnal
        size of a pixel on sky. default: 0.00125
        note: fieldOfView/pixelScale must be the number of pixels accross the pupil file(400)
    saveIntermediateComplexWavefront: boulean. optionnal
        Do you want to save each intermediate complex wavefront, default : False
    pathFileComplex: string, optionnal
        Where to save the complexe wavefront
    primaryFits: string, optionnal
        path to the primary segment map
    primaryPistRms: float, optionnal
        rms piston error on the primary mirror in meters. 
    
    
    Returns
    -------
    imageFileName : str
        full name of the image generated at the end of this function.
    '''

    image=np.zeros((nPix,nPix))
    rng = np.random.default_rng(seed)
    amp= rng.normal(0, rmsIslandEffect, 6) #mm
    if rmsIslandEffect>0:
        amp *= rmsIslandEffect/np.std(amp)
        amp -= np.mean(amp)
    staticPist = Mod_segment(pupilMap,amp)
    # print('standard deviation: {}'.format(np.std(staticPist[pupilMap>0])))
    
    # plt.figure()
    # plt.title('static pist')
    # plt.imshow(staticPist)
    # plt.colorbar()
    # plt.axis('off')
    
    if primaryFits is not None:
        prim = ELTprimary(primaryFits)
    
    if esoData:
        if esoPsfTipTilt != None:
            ttData = fits.open(path2esoData+esoPsfTipTilt)
            rng = np.random.default_rng(seed)
            start=rng.integers(0,ttData[0].data.shape[1]-nScreen-1)
            
            ttSequence = ttData[0].data[:,start:start+nScreen]
            for tt in range(2):
                ttSequence[tt] *= rmsWindshake/np.std(ttSequence[tt])
                ttSequence[tt] -= np.mean(ttSequence[tt])
            
            ttData.close()
        if nScreen<2:
            print('warning: not enough screens to rescale the tip tilt, omiting tiptilt')
            ttSequence = np.array([[0],[0]])
    
    # for i in np.arange(0,nScreen):
    for i in range(nScreen):
        print('\r{:04}'.format(i),end = '')
       
        
        atmScreen = po.FITSOpticalElement( transmission=pupilFileName
                                       ,opd=atmFilePrefix+'{:.0f}.fits'
                                       .format(i+300),opdunits='meters'
                                      ,planetype=po.poppy_core.PlaneType.pupil)
        
        # tmp = atmScreen.opd.copy()
        
        #scale the atmosphere if a different 
        atmScreen.opd*=atmWorseningFactor
        
        #remove AO residual piston opd
        for s in range(6):
            atmScreen.opd[pupilMap==s+1]-=np.mean(atmScreen.opd[pupilMap==s+1])
        # plt.figure()
        # plt.title('residue')
        # plt.imshow(atmScreen.opd)
        # plt.colorbar()
        # plt.axis('off')
        
        #create the dynamic piston
        jit = rng.normal(0,rmsSegmentJitter,6)
        if rmsSegmentJitter>0:
            # jit *= rmsIslandEffect/np.std(jit)
            jit -= np.mean(jit)
        jitterPist = Mod_segment(pupilMap, jit)
        
        # plt.figure()
        # plt.title('dynamicPost')
        # plt.imshow(jitterPist)
        # plt.colorbar()
        # plt.axis('off')
        
        if primaryFits is not None:
            pist = rng.normal(0,1,len(prim.segmentlist))*primaryPistRms
            for s in range(len(prim.segmentlist)):
                prim.set_segment(prim.segmentlist[s],pist[s])
            wave = po.Wavefront(npix = 400)
            primopd = prim.get_opd(wave)
        else:
            primopd = np.zeros((atmScreen.opd.shape))
        
        # plt.figure(1)
        # plt.clf()
        # plt.imshow(primopd)
        
        #add the contribution of the static and dynamic piston
        atmScreen.opd+=staticPist+jitterPist+primopd
        
        # plt.figure(2)
        # plt.clf()
        # plt.imshow(atmScreen.opd)
        
        #setup the windshake
        if esoData and esoPsfTipTilt!=None:
            tip = ttSequence[0,i]
            tilt = ttSequence[1,i]
        else:
            rWS = rng.normal(0,rmsWindshake)
            dirWS = rng.uniform(-np.pi,np.pi)
            
            tip = rWS*np.cos(dirWS)
            tilt = rWS*np.sin(dirWS)
        
        wfe = po.ZernikeWFE(radius=39/2,
                            coefficients=[0,tip,tilt ],
                            aperture_stop=True)
        
        #set up the optica system for calculation
        osys = po.OpticalSystem(oversample = 4)
        osys.add_pupil(atmScreen)
        osys.add_pupil(wfe)
        osys.add_detector(pixelscale=pixelScale, fov_arcsec=fieldOfView)
        
        #calculation
        # psf = osys.calc_psf(wavelength, display_intermediates=True)#for debugging
        if saveIntermediateComplexWavefront:
            psf,compWf = osys.calc_psf(wavelength, return_intermediates = True)
            
            hdu =fits.PrimaryHDU(np.array([np.real(compWf[-1].wavefront)
                                           ,np.imag(compWf[-1].wavefront)]))
            
            hdul = fits.HDUList(hdu)
            hdul[0].header['ATMFACTO'] = (atmWorseningFactor ,'factor applied to the atmospheric screen')
            hdul[0].header['ATMSERIE'] = atmFilePrefix
            hdul[0].header['RNGSEED'] = seed
            hdul[0].header['RMSSTAT'] = rmsIslandEffect
            hdul[0].header['RMSSEGJI'] = rmsSegmentJitter
            hdul[0].header['RMSWS'] = rmsWindshake
            hdul[0].header['MAGNITUD'] = outputFileSufix
            hdul[0].header['SCREENNB'] = i+300
            if not(os.path.exists(pathFileComplex)):
                os.mkdir(pathFileComplex)
            hdul.writeto(pathFileComplex+outputFileSufix+'screenNb{}.fits'.format(i+300))
            del hdul[0].data
            
        else:
            psf = osys.calc_psf(wavelength)
        
        # plt.figure()
        # plt.imshow(np.log(psf[0].data),cmap = 'gist_heat')
        # plt.axis('off')
        image += psf[0].data
    
    print('\nfinished')
    
    # plt.figure()
    # plt.imshow(np.log(image),cmap = 'gist_heat')
    # plt.axis('off')

    psf[0].data = image.copy()
    psf[0].header['ATMFACTO'] = (atmWorseningFactor ,'factor applied to the atmospheric screen')
    psf[0].header['ATMSERIE'] = atmFilePrefix
    psf[0].header['RNGSEED'] = seed
    psf[0].header['RMSSTAT'] = rmsIslandEffect
    psf[0].header['RMSSEGJI'] = rmsSegmentJitter
    psf[0].header['RMSWS'] = rmsWindshake
    psf[0].header['MAGNITUD'] = outputFileSufix
    filename1 = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")+'_'+outputFileSufix+'_{:.0f}nm.fits'.format(wavelength*10**9)
    psf.writeto(filename1)
    
    print('sucessfull write')

    return filename1


if __name__ == '__main__':
    ppMap = buildPupil()
    
    # outFile = SCAOSim(1000*10**-9, 1600, 1000, 00*10**-9, 00*10**-9, 100*10**-9, 0.0
    #                   , 12345, ppMap,atmFilePrefix = 'screen_I=8/screenI=8'
    #                   ,outputFileSufix='test')
    
    outFile = SCAOSim(1000*10**-9, 1600, 1000, 00*10**-9, 00*10**-9, 000*10**-9, 0.0
                      , 12345, ppMap,atmFilePrefix = 'SEEING/seeing1/ResidualOPD_seeing_1_iteration_'
                      ,outputFileSufix='test',esoData=True,
                      path2esoData='C:/Users/anne.cheffot/simulations/python/ANDES-SCAO/',
                      esoPsfTipTilt='Time_hist_TT_wind10ms_TelZ45M2_290sec_500Hz.fits',
                      primaryFits='primaryMap.fits',primaryPistRms=100*10**-9)




