MR_libs
=======
### Description

These are MATLAB libraries for post-processing, analysing and simulating Magnetic Resonance Spectroscopy and Imaging (MRS &amp; MRI) data.

### Configuration

Data files & header files are required to be added to MATLAB path.  

### File List
**MRS_lib:** contains functions for post-processing, analysing and simulating Magnetic Resonance Spectroscopy data
* io/ 
  * **mrs_readSPAR.m**  : reads .SPAR Philips MRS header file
  * **mrs_readSDAT.m**  : reads .SDAT Philips MRS data file
  * **mrs_readLIST.m**  : reads .list Philips MRS header file produced by delayed reconstruction
  * **mrs_readDATA.m**  : reads .data Philips MRS data file produced by delayed reconstruction
  * **mrs_readSIN.m**   : reads .data Philips raw MRS data file
  * **mrs_readRAW.m**   : reads .data Philips raw MRS data file
  * **mrs_readLcmodelRAW.m** : reads LCModel input .RAW file which contains time domain data of each metabolite spectrum
  * **mrs_readLcmodelBASIS.m**   : reads the LCModel input .basis file which contains the basis set of model metabolite spectra 
  * **mrs_readLcmodelTABLE.m**   : reads the metabolite absolute and relative concentration and their SDs from the LCModel output .table file 
  * **mrs_readLcmodelCOORD.m**   : reads the LCmodel output .coord file, which contains the coordinates of all curves on the one-page output
  * **mrs_readJmruiTXT.m**   : reads .txt MRS data file from jMRUI
  * **mrs_writeSDAT.m** : writes MRS data to Philips .SDAT file
  * **mrs_editSPAR.m**  : edits Philips .SPAR file
* postprocess/
  * **mrs_truncate.m**  : truncates points from the end of spectra
  * **mrs_zerofill.m**  : fills zeros to the end of spectra
  * **mrs_apod.m**      : applies line-broadening filter 
  * **mrs_fft.m**       : applies fourier transformation 
  * **mrs_ifft.m**      : applies inverse fourier transformation 
  * **mrs_rephase.m**   : rephases spectra with specified phase value  
  * **mrs_zeroPHC.m**   : applies automatic zero-order phase correction to a spectrum
  * **mrs_manualzeroPHC.m**   : allows users to manually apply zero-order phase correction of a spectrum
  * **mrs_firstPHC.m**  : applies automatic first-order phase correction to a spectrum
  * **mrs_realign.m**   : aligns the target peaks in the spectra
  * **mrs_average.m**   : calculates averaged data
  * **mrs_lowpass.m**   : filters out high frequency components above specified frequency
  * **mrs_highpass.m**  : filters out low frequency components below specified frequency     
* simulation/  
  * **mrs_PRESS.m**      : simulates signal acquired using Position Resolved Spectroscopy (PRESS)
  * **mrs_STEAM.m**      : simulates signal acquired using Stimulated Echo Acquisition Mode (STEAM) 
  
  ![Alt text](https://raw.github.com/chenkonturek/MR_libs/master/Images/pulse_seq.jpg)  
  * **mrs_sLASER.m**     : simulates the signal acquired using semi-localised by adiabatic selective refocusing(sLASER) sequence 
  * **mrs_ISISscheme.m** : demonstrates how Image Selective in vivo Spectroscopy (ISIS) works
  * **mrs_simulateFID.m**: simulates a Free Induction Decay (FID) or Half-Echo. 
* utils/
  * **mrs_selectPeakrange.m** : allows users to manually pick the peak range interactively
  * **mrs_findPeak.m**        : locates the highest positive peak or lowest negative peak
  * **mrs_fitPeak.m**         : fits a peak in the given range of a spectrum with a lorenztian curve  
  * **mrs_lorentzFit.m**      : fits data with a lorenztian function by minimising the squared error
  * **mrs_lorentzFun.m**      : defines the lorentzian function
  * **mrs_points2Hz.m**       : converts unit from points to Hz
  * **mrs_points2ppm.m**      : converts unit from points to ppm
  * **mrs_ppm2Hz.m**          : converts unit from ppm to Hz
  * **mrs_plotSpectra.m**     : displays spectra 
  * **mrs_viewCSI.m**         : displays a spectrum from a selected voxel      
  * **mrs_rot90.m**           : rotates the spectra images 90 degree clockwise
  * **mrs_calTemp.m**         : calculates the temperature based on chemical shift difference of the water resonance and the temperature-independent reference resonance

**MRI_lib** contains functions for post-processing, analysing and simulating Magnetic Resonance Imaging data
* io/
  * **mri_readHDR.m**  : reads .hdr MRI ANALYZE 7.5 header file 
  * **mri_readIMG.m**  : reads .img MRI ANALYZE 7.5 data file 
  * **mri_writeIMG.m** : writes data to .img file 
  * **mri_readPAR.m**  : reads .PAR Philips MRI header file
  * **mri_readREC.m**  : reads .REC Philips MRI data file
* utils/
  * **mri_IRcurve.m**        : calculates the inversion recovery (IR) curve

  ![Alt text](https://raw.github.com/chenkonturek/MR_libs/master/Images/IR.png)  

  * **mri_absIRcurveFit.m**  : fits the data to absolute inversion recovery curve 
  * **mri_locateSVOI.m**     : creates a mask for locating the spectroscopic VOI in the axial MR images
  * **mri_dispSVOI.m**       : displays the spectroscopic VOI on top of the axial MR images

**NMR_lib:** contains functions for NMR simulation 
  * **nmr_bloch.m** : defines full Bloch equations in rotating frame
  * **nmr_calT2star.m** : calculates T2*
  * **nmr_getGamma** : returns gyromagnetic ratio value for different nuclei

**Examples:** contains example scripts. (Please email me if you want the data files) 
  * **example1.m** : demonstrates how to use MR_libs to post-process MRS data.
  * **example2.m** : demonstrates how to use MR_libs to locate spectroscopic voxel in MR images. 
  * **example3.m** : demonstrates how to use MR_libs to to simulate tye half-echo signal acquired using PRESS and STEAM sequences.respectively.
  * **example4.m** : demonstrates how to use MR_libs to do Bloch Equation simulation. 
  * **example5.m** : demonstrates how to use MR_libs to simulate magnetisation profile produced by an RF pulse.
  * **example6.m** : demonstrates how to use MR_libs to estimate T1.
  * **example7.m** : demonstrates how to use MR_libs to displays a spectrum from a selected voxel.

### Acknowledgements

I would like to thank Professor Penny Gowland and Dr. Susan Francis for their supervision. 
I would also like to express my enormous thanks to Dr. Mary Stephensons and Dr. Emma Halls for their help and contributions.      


