"""
Created on Wed Apr  1 14:28:29 2020

@author: esthersoriahernandez

This function calculate the center of gravity of the input image, discarding the pixels
that are below a given threshold
Inputs:
	*Image or area 
	*Thershold (absolut value)
Outputs:
	*COG coordinate x
	*COG coordinate y

"""
def cog(input, threshold):
    import numpy as np 
    #
    
    height = np.size(input,0)
    width = np.size(input,1)
    
    sum_numx=0
    sum_numy=0
    sum_den=0
    pixel=np.zeros((height,width))
    for i in range(1,width+1):
        for j in range(1,height+1):
            pixel=input[j-1,i-1]
            if (pixel>threshold):
                pixel=pixel-threshold
                sum_numx=sum_numx+pixel*(i-1)
                sum_numy=sum_numy+pixel*(j-1)
                sum_den=sum_den+pixel

    cx=sum_numx/sum_den
    cy=sum_numy/sum_den
    return cx,cy
