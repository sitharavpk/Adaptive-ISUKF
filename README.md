# Adaptive Importance Sampling Unscented Kalman Filter based SAR Image Super Resolution

This is a Matlab implementation of SAR Image Super Resolution. SAR images being inherently affected by speckle 
noise fails on using natural image super-resolution methods. The code presented here simultaneously denoises 
and super-resolves to a single high-resolution frame from multiple low-resolution images.

Authors: Sithara Kanakaraj, Madhu S. Nair and Saidalavi Kalady

This code is the Matlab implementation of the paper

Sithara Kanakaraj, Madhu S. Nair and Saidalavi Kalady, “Adaptive Importance Sampling Unscented Kalman Filter based SAR Image Super Resolution”, Computers and Geosciences, Elsevier, Vol. 133, Article No. 104310, December 2019. 
```
DOI: https://doi.org/10.1016/j.cageo.2019.104310
```

To run the code, please use the main function Adaptive_ISUKF.m.

## Input:
    The test images are present in the data folder. Each folder consists of 
    synthetic and real SAR images. Each input LR image consists 16 low-resolution images. 
    
   Format:
         - filePath: the file name of the input 
            Example:    filePath = 'data\synthetic\synthetic1.mat';

## Output:
    A high-resolution image magnified to a factor of 2 is presented as the output along 
    with the values of quality assessment metrics.

        - HR: the super-resolved image

## External codes

    1. For Image Registration: Manuel Guizar-Sicairos, Samuel T. Thurman, and James R. Fienup, 
    "Efficient subpixel image registration algorithms," Opt. Lett. 33, 156-158 (2008).
    
    2. For Noise estimation: Xinhao Liu, Masayuki Tanaka and Masatoshi Okutomi, 
    "Single-Image Noise Level Estimation for Blind Denoising Noisy Image", IEEE Transactions 
    on Image Processing, Vol.22, No.12, pp.5226-5237, December, 2013.

    3. For Structural Similarity Index Metric (SSIM): Z. Wang, A. C. Bovik, H. R. Sheikh, and 
    E. P. Simoncelli, "Image quality assessment: From error visibility to structural similarity,"
    IEEE Transactios on Image Processing, vol. 13, no. 4, pp. 600-612, Apr. 2004.

    4 For Feature Similarity Index Metric (FSIM): Lin Zhang, Lei Zhang, Xuanqin Mou, and David Zhang,
    "FSIM: a feature similarity index for image qualtiy assessment", IEEE Transactions on Image 
    Processing, vol. 20, no. 8, pp. 2378-2386, 2011.

## Citation
Please cite our paper if you find the software useful for your work.

```
@article{kanakaraj2019adaptive,
  title={Adaptive Importance Sampling Unscented Kalman Filter based SAR image super resolution},
  author={Kanakaraj, Sithara and Nair, Madhu S and Kalady, Saidalavi},
  journal={Computers \& Geosciences},
  volume={133},
  pages={104310},
  year={2019},
  publisher={Elsevier}
}
```    
Disclaimer: The assessment metric values in the paper are the best results for the ideal cases. It may change at every execution.

