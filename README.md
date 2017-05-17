# VAMP Deconvolution

Original release date : 12/5/17

Reference #1          : "Vector Approximate Message Passing"
Authors               : Philip Schniter, Sundeep Rangan, Alyson K. Fletcher
Download              : https://arxiv.org/abs/1610.03082

Reference #2          : "Generalized Expectation Consistent Signal Recovery for Nonlinear Measurements"
Authors               : Hengtao He, Chao-Kai Wen, Shi Jin
Download              : https://arxiv.org/abs/1701.04301

Contrubutors : Alexis Gladiline Bozio, Florent Krzakala, Anthony Maggs, Andre Manoel, Eric Tramel 

Contents
---------------------------------------------------------------------------
scripts: 


    CS_Demo_Vamp_on_convoluted_images.m : Deconvolution of a blurred and downsampled image with the VAMP algorithm with the use of the
    selective inversion method and BM3D denoiser    

    Demo_compare_timing_vamp.m : Deconvolution of a blurred and downsampled toy image (only 0 and 1) with the Vector Approximate 
    Message Passing algorithm. Comparison in timing between the selected inversion and the use of Woosbury's identity
    
    
    Demo_streaming_vamp.m : Deconvolution of a blurred and downsampled toy image (only 0 and 1) with the Vector Approximate Message 
    Passing algorithm using the selected inversion method and by dividing the image by block
    
    Demo_Super_resolution_VAMP.m : Super_resolution reconstruction with Vamp and the Generalized Expectation Consistent Signal Recovery 
    for Nonlinear Measurements.


This repository includes the BM3D package.
The latest versions of this package can be found at:
    BM3D: http://www.cs.tut.fi/~foi/GCF-BM3D/

This repository is calling the Selinv implementation of the Selected Inversion algorithm, which computes the selected elements of the inverse of a general sparse symmetric matrix. 
The latest versions of this package can be found at:
    https://math.berkeley.edu/~linlin/software.html
