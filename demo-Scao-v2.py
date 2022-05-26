"""
Created on Fri Jan 28 11:45:35 2022

@author: anne.cheffot
@author: EstherSoria



(ALC)okay here is my suggestion This script needs to contain only PSF computation
No radial profile drawing no display of PSF (at the end, of course for debugging
purposes you can do what ever). The function would have the following prototype

def SCAOSim(wavelength,nPix,nScreen,rmsIslandEffect,rmsSegmentJitter,rmsWindshake
            ,seed,atmWorseningFactor,pupilFileName='Tel-Pupil.fits', atmFilePrefix='screen'
            , outputFileSufix='I=8-2200nm.fits' ):
    '''
    Parameters
    ----------
    wavelength : float
        wavelength of the final image in meters
    nPix : int
        number of pixels across one image
    nScreen : int
        number of time light is propagated through a aberated screen.
    rmsIslandEffect : float
        root mean square error of the island effect in meters
    rmsSegmentJitter : float
        root mean square of the segment jitter in meters
    rmsWindshake : float
        root mean square of the tip tilt induced by the winshake in meters
    seed : int
        seed for the randome number generator
    atmWorseningFactor : float
        a number to change the atmospheric residuals amplitude (no unit)
    pupilFileName : string, optional
        The name of the file in which the pupil shape is stored. The default is 
        'Tel-Pupil.fits'.
    atmFilePrefix : string, optional
        The atmosheric screens should have a prefix, for example 'screen' and be 
        followed by an iterger. The default is 'screen'.
    outputFileSufix : string, optional
        The output file should be date_time followed by a suffix.
        The default is 'I=8-2200nm.fits'.

    Returns
    -------
    It could return the name of the file saved at the end

    '''


"""
import numpy as np
import matplotlib.pyplot as plt
import poppy as po
# from astropy.io import fits
# from astropy import units as u
import datetime
from pupil_ANDES import *
from cog import *
# from rad_distrub import *#This library does not exist for me is it one of yours?

nPix = 400

wavelength = 2.2*10**-6
# startScreen = 200
# worseningFactor = 0.0#diffraction limited
worseningFactor = 1.0 #I = 8
# worseningFactor = 1.6 #I = 11.5
# worseningFactor = 2.16 #I = 15
# worseningFactor = 2.57 #I = 
nScreen = 200

rmsIslandEffect=500*10**-9
rmsSegmentJitter = 10*10**-9
rmsWindshake = 100*10**-9
pupilFileName='Tel-Pupil.fits'
atmFilePrefix='screen'
outputFileSufix = 'I=8-2200nm.fits'

#here I would like to suggest a simplification: you do not need to specify 
#segm and amp you only need amp that must be length 6 and contain the piston of 
#each segments (even when that segment has no phasing error). 

segm=[1,2,3,4,5,6]
#(ALC) here we need to introduce 2 things:repeatability and controlability
#we want to be able to repeat configuration that worked/did not worked we do this
#through seeds. here is an example code of how to do this.
# seed = 12345#if you set a seed you will always get the same numbers out of the 
# #random number generator, change the seed to get a different set
# rmsIslandEffect=500*10**-9#this introduces controlability you know what rms you will get 
# #out of your piston generation
# rng = np.random.default_rng(seed)
# amp = rng.normal(0,rms,6) # this create random numbers centered on 0 with a 
#sigma of rms and it creates 6 of them

amp= np.random.randint(-2, 2, 6, dtype=int) #um
rms = np.sqrt(np.mean(amp**2))
print(rms) #the rms generated in um
#(ALC) I guess misalignment means cophasing error for you?
ampli=Mod_segment(amp,segm,1) #generate the opd fits file for the misaligment
#(ALC) we need to talk about flags! This is very bad practice
alig=0 #1 if you want to add the contributions of the vibrations
vibra=0 #1 if you want to add the contributions of the misalignment in the segments
DL=0 #1 if you want to compare with the difraction limit results



image = np.zeros((nPix,nPix))
for i in range(nScreen):
    print('\r{:04}'.format(i),end = '')
    #(ALC)The two following if are not needed if you want the DL you put the 
    # worseningFactor to 0
    if DL!=1:
        atmScreen = po.FITSOpticalElement( transmission=pupilFileName
                                       ,opd=atmFilePrefix+'{:.0f}.fits'
                                       .format(startScreen+i),opdunits='meters'
                                      ,planetype=po.poppy_core.PlaneType.pupil)
    if DL==1:
        atmScreen = po.FITSOpticalElement( transmission=pupilFileName
                                       ,opd='Zeros.fits'
                                      ,planetype=po.poppy_core.PlaneType.pupil)
    
    atmScreen.opd*=worseningFactor
    piston_segm=po.FITSOpticalElement(transmission=pupilFileName, opd='Segm_pupil.fits'
                                      ,opdunits='nm',planetype=po.poppy_core.PlaneType.pupil)
    osys = po.OpticalSystem(npix = nPix,oversample = 4)
    osys.add_pupil(atmScreen) #adding the residuals
    
    
    
    amp= np.random.normal(0, rmsSegmentJitter, 6) #um
    #(ALC)for repeatability use rng.normal
    ampli=Mod_segment(amp,segm,2)
    
    
    #wfe = po.ZernikeWFE(radius=15*u.m,#(ALC)this seems too small, why 15? and not 38.542//2)
                        #coefficients=[0,rng.normal(0,rmsWindshake), rng.normal(0,rmsWindshake)],
                        #aperture_stop=True)
  
    #osys.add_pupil(wfe) simlating a tiptil move for the windshake
    
    if alig==1:    #adding the misaligment
        osys.add_pupil(piston_segm)
    if vibra==1:   #adding the vibration
        piston_vibra=po.FITSOpticalElement(transmission=pupilFileName, opd='Vibra_pupil.fits'
                                      ,opdunits='nm',planetype=po.poppy_core.PlaneType.pupil)
        osys.add_pupil(piston_vibra)
    osys.add_detector(pixelscale=0.005, fov_arcsec=0.5)
    # psf = osys.calc_psf(1.6e-6, display_intermediates=True)#for debugging
    psf = osys.calc_psf(wavelength)
    image += psf[0].data #appending the data to simulate a high exposition

print('\nfinished')

#Save the results with the folliwing structure: date+mag+wavelength
psf[0].data = image.copy()
psf[0].header['ATMFACTO'] = (worseningFactor,'factor applied to the atmospheric screen')
filename1 = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
psf.writeto(filename1+outputFileSufix)
print('sucessfull write')
# plt.figure()
# plt.imshow(np.log(image))
#plt.figure()
#wfe.display(what='both');
#Analysis PSF, save in different variables the results for each situation
#(ALC)I strongly advise separating the radial profile drawing in a separated 
#script. I'll clean up my dedicated script so that it is 
#this does not work for me because I cannot find rad_distrub
# if DL==1:
#     imageDL=image
#     centx,centy=cog(imageDL,1e-06)
#     radDL,pixDL=RadialProfile(imageDL,centy, centx)
#     norm_fact=np.max(radDL)
#     radDL=radDL/norm_fact
# else:
#     if vibra!=1 and alig!=1:
#         centx,centy=cog(image,1e-06)
#         rad,pix=RadialProfile(image,centy, centx)
#         norm_fact=np.max(rad)
#         rad=rad/norm_fact
#     if vibra!=1 and alig==1:
#         image1=image
#         centx,centy=cog(image1,1e-06)
#         rad1,pix1=RadialProfile(image1,centy, centx)
#         norm_fact=np.max(rad1)
#         rad1=rad1/norm_fact
#     if vibra==1 and alig!=1:
#         image3=image
#         centx,centy=cog(image3,1e-06)
#         rad3,pix3=RadialProfile(image3,centy, centx)
#         norm_fact=np.max(rad3)
#         rad3=rad3/norm_fact
#     if vibra==1 and alig==1:
#         image2=image
#         centx,centy=cog(image2,1e-06)
#         rad2,pix2=RadialProfile(image2,centy, centx)
#         norm_fact=np.max(rad2)
#         rad2=rad2/norm_fact

#(ALC) in python it is bad practice to put commented block in three quote blocks
"""
#plotingt the results:
plt.figure()
plt.yscale("log")
plt.title('Radial distribution, I=8, WL 2,2um', fontsize=15)
plt.xlabel('Distance from PSF center (mas)', fontsize=12)
plt.ylabel('Frequency', fontsize=12)
plt.plot(pixDL[0:80]*500/400, radDL[0:80],  label='psf DL')
plt.plot(pix[0:80]*500/400, rad[0:80],  label='psf + res')
#plt.plot(pix1[0:80]*500/400, rad1[0:80], label='psf + pist + res')
#plt.plot(pix3[0:80]*500/400, rad3[0:80], label='psf + vibra + res'),
#plt.plot(pix2[0:80]*500/400, rad2[0:80], label='psf pist + vibra + res')
plt.grid()
plt.legend(loc=1)
plt.show()



#ts=aotools.azimuthal_average(np.abs(image))

"""
