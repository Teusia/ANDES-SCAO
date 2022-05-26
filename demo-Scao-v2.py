"""
Created on Fri Jan 28 11:45:35 2022

@author: anne.cheffot
@author: EstherSoria
"""
import numpy as np
import matplotlib.pyplot as plt
import poppy as po
from astropy.io import fits
from astropy import units as u
import datetime
from pupil_ANDES import *
from cog import *
from rad_distrub import *

#if you want primary mirror segmentation use this. but not recomended
# seglist = [a for a in range((3*4**2)+(3*4)+1,(3*15**2)+(3*15)+1)]
# seglist.extend([a for a in range(723,736)])
# seglist.extend([a for a in range(739,752)])
# seglist.extend([a for a in range(755,768)])
# seglist.extend([a for a in range(771,784)])
# seglist.extend([a for a in range(787,800)])
# seglist.extend([a for a in range(803,816)])

# seglist.extend([a for a in range(821,831)])
# seglist.extend([a for a in range(838,848)])
# seglist.extend([a for a in range(855,865)])
# seglist.extend([a for a in range(872,882)])
# seglist.extend([a for a in range(889,899)])
# seglist.extend([a for a in range(906,916)])

# ap = po.MultiHexagonAperture(rings = 15, flattoflat = 1.239*u.m,gap=0.004*u.m,
#                               segmentlist = seglist,
#                               rotation = 0)
# pupilStop = po.CircularAperture(radius=38.542/2)
# sec = po.SecondaryObscuration(secondary_radius= 6.5*u.m, n_supports = 6,
#                               support_width = 0.54*u.m,rotation = 30)
# pupil = po.CompoundAnalyticOptic(opticslist=[
#     ap,pupilStop,
#     sec
#     ],name = 'ELT pupil')

image = np.zeros((400,400))
wavelength = 2.2*10**-6
startScreen = 1
# worseningFactor = 0.0#diffraction limited
worseningFactor = 1.0 #I = 8
# worseningFactor = 1.6 #I = 11.5
# worseningFactor = 2.16 #I = 15
# worseningFactor = 2.57 #I = 14
segm=[1,2,3,4,5,6]
amp= np.random.randint(-2, 2, 6, dtype=int) #um
rms = np.sqrt(np.mean(amp**2))
print(rms) #the rms generated in um
ampli=Mod_segment(amp,segm,1) #generate the opd fits file for the misaligment
alig=0 #1 if you want to add the contributions of the vibrations
vibra=0 #1 if you want to add the contributions of the misalignment in the segments
DL=0 #1 if you want to compare with the difraction limit results

for i in range(200):
    print('\r{:04}'.format(i),end = '')
    if DL!=1:
        atmScreen = po.FITSOpticalElement( transmission='Tel-Pupil.fits'
                                       ,opd='part3-screen{:.0f}.fits'
                                       .format(startScreen+i),opdunits='meters'
                                      ,planetype=po.poppy_core.PlaneType.pupil)
    if DL==1:
        atmScreen = po.FITSOpticalElement( transmission='Tel-Pupil.fits'
                                       ,opd='Zeros.fits'
                                       .format(startScreen+i)
                                      ,planetype=po.poppy_core.PlaneType.pupil)
    
    tmp = atmScreen.opd.copy()
    
    atmScreen.opd*=worseningFactor
    piston_segm=po.FITSOpticalElement(transmission='Tel-Pupil.fits', opd='Segm_pupil.fits'.format(startScreen+i)
                                      ,opdunits='nm',planetype=po.poppy_core.PlaneType.pupil)
    osys = po.OpticalSystem(npix = 400,oversample = 4)
    osys.add_pupil(atmScreen) #adding the residuals
    amp= np.random.normal(0, 0.01, 6) #um
    ampli=Mod_segment(amp,segm,2)
   
    #wfe = po.ZernikeWFE(radius=15*u.m,
                        #coefficients=[np.random.random()*1e-7, np.random.random()*1e-7, 1e-7, 0e-6, 0, 0e-6, 0e-6],
                        #aperture_stop=True)
  
    #osys.add_pupil(wfe) simlating a tiptil move for the windshake
    
    if alig==1:    #adding the misaligment
        osys.add_pupil(piston_segm)
    if vibra==1:   #adding the vibration
        piston_vibra=po.FITSOpticalElement(transmission='Tel-Pupil.fits', opd='Vibra_pupil.fits'.format(startScreen+i)
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
if worseningFactor==0.0 :
    psf.writeto(filename1+'I=8-{:.0f}nm-DL.fits'.format(wavelength*10**9),overwrite = True)
elif worseningFactor==1.0 :
    psf.writeto(filename1+'-I=8-{:.0f}nm.fits'.format(wavelength*10**9),overwrite = True)
elif worseningFactor == 1.6: 
    psf.writeto(filename1+'-I=11_5-{:.0f}nm.fits'.format(wavelength*10**9),overwrite = True)
elif worseningFactor == 2.16:
    psf.writeto(filename1+'-I=15-{:.0f}nm.fits'.format(wavelength*10**9),overwrite = True)
elif worseningFactor == 2.57:
    psf.writeto(filename1+'-I=14-{:.0f}nm.fits'.format(wavelength*10**9),overwrite = True)
print('sucessfull write')
plt.figure()
plt.imshow(np.log(image))
#plt.figure()
#wfe.display(what='both');
#Analysis PSF, save in different variables the results for each situation
if DL==1:
    imageDL=image
    centx,centy=cog(imageDL,1e-06)
    radDL,pixDL=RadialProfile(imageDL,centy, centx)
else:
    if vibra!=1 and alig!=1:
        centx,centy=cog(image,1e-06)
        rad,pix=RadialProfile(image,centy, centx)
    if vibra!=1 and alig==1:
        image1=image
        centx,centy=cog(image1,1e-06)
        rad1,pix1=RadialProfile(image1,centy, centx)
    
    if vibra==1 and alig!=1:
        image3=image
        centx,centy=cog(image3,1e-06)
        rad3,pix3=RadialProfile(image3,centy, centx)
    
    if vibra==1 and alig==1:
        image2=image
        centx,centy=cog(image2,1e-06)
        rad2,pix2=RadialProfile(image2,centy, centx)


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
