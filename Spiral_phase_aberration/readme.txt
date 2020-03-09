--------------------------------------------------------------------------------------------------------------------------

 Code demo for "Shack-Hartmann Fourier ptychography microscopy (SHFPM)"
 Public release v1.0 (Mar 9th, 2020) 

------------------------------------------------------------------------------------------------------------------------------------
Contents in EPRY folder 
this is a conventional EPRY code, which can be find in Guoan Zheng, Fourier Ptychographic Imaging A MATLAB® tutorial, Morgan & Claypool Publishers
------------------------------------------------------------------------------------------------------------------------------------
***************** main functions *********************
*) ini_enviroment.m           : This function simulates the experimental environment;
*) EPRY_getimage.m            : This function simulates captured low resolution images of FPM using microlens arrays;
*) EPRY_reconstruction.m      : This function runs the conventional  EPRY algorithm;


------------------------------------------------------------------------------------------------------------------------------------
Contents in SHFPM folder
------------------------------------------------------------------------------------------------------------------------------------
***************** main functions *********************
*) Demo.m                      : This demo does the simulation of the Shack-Hartmann Fourier ptychographic microscopy technique in the case with spiral phase aberrations
*) ini_enviroment.m            : This function simulates the experimental environment;
*) shFPM_getSubimage.m         : This function simulates captured low resolution images of FPM using microlens arrays;
*) shFPM_reconstruction.m      : This function runs the proposed SHFPM algorithm;

***************** other functions ********************
*) arrayfft2.m   : This function performs array-like Fourier tranformation corresponding to the function of microlens array


If you have any comment, suggestion, or question, please feel free to contact Shuhe Zhang at shuhe.zhang@maastrichtuniversity.nl
