# Super-Resolution-Imaging-with-Message-Passing-Algorithms

Original release date : 12/5/17

Reference #1          : "Vector Approximate Message Passing"
Authors               : Philip Schniter, Sundeep Rangan, Alyson K. Fletcher
Download              : https://arxiv.org/abs/1610.03082



Contents
---------------------------------------------------------------------------
scripts:
    CDemo.m: Recover a compressively sampled 1D signal with Haar wavelet sparsity based (V)AMP and NLM-(V)AMP.
    

functions:
    AMP.m: Reconstructs sparse, compressively sampled signals using AMP.
    VAMP.m: Reconstructs sparse, compressively sampled signals using VAMP.
    DAMP.m: Performs D-AMP reconstruction of a compressively sampled signal. The string "denoiser" selects which denoiser to use.
    DVAMP.m: Performs D-VAMP reconstruction of a compressively sampled signal. The string "denoiser" selects which denoiser to use.
	DprGAMP: Performs D-prGAMP compressive phase retrieval of a signal. The string "denoiser" selects which denoiser to use.
    DIT.m: Performs D-IT reconstruction of a compressively sampled signal. The string "denoiser" selects which denoiser to use.
    DAMP_oneIter.m: Performs a single iteration of D-AMP.  Used to generate state evolution and qqplots.
    DIT_oneIter.m: Performs a single iteration of D-IT.  Used to generate state evolution and qqplots.
    DAMP_SE_Prediction.m: Computes the next predicted state evolution of a D-AMP algorithm.
    denoise.m: Denoises an incoming signal.  Which denoiser is used depends on the string "denoiser".  Currently supports Gaussian filtering, bilateral filtering, NLM, BLS-GSM, BM3D, and BM3D-SAPCA.  Add your own denoising algorithms here.



Packages
---------------------------------------------------------------------------
This download includes the BM3D, BLS-GSM, NLM, and Rice Wavelet Toolbox packages.
The latest versions of these packages can be found at:
    BM3D: http://www.cs.tut.fi/~foi/GCF-BM3D/
    BLS-GSM: http://decsai.ugr.es/~javier/denoise/software/index.htm
    NLM: http://www.mathworks.com/matlabcentral/fileexchange/27395-fast-non-local-means-1d--2d-color-and-3d
    Rice Wavelet Toolbox (RWT): https://github.com/ricedsp/rwt
