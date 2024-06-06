% Motion Correction for fMRI phase data
% V Malekian, FIL Physics, March 2023

%% Prepare magnitude and phase data (merging 3D nifti files into 4D format) 
directory_m='/export/home/vmalekian/mag/'; %% change the name of directory to your local one
unix(['fslmerge -t ' directory_m 'mg ' directory_m 'f*.nii' ])
directory_p='/export/home/vmalekian/phase/'; %% change the name of directory to your local one
unix(['fslmerge -t ' directory_p 'ph ' directory_p 'f*.nii' ])

%% Rescale phase data to 0 to 2pi range
unix(['fslmaths ' directory_p 'ph -mul 3.14159 -div 2048 ' directory_p 'ph_rad -odt float'])

%% Extract real and imaginary data from magnitude and phase 
unix(['fslcomplex -complexpolar ' directory_m 'mg ' directory_p 'ph_rad ' directory_m 'cx 0 168']) % In this example fMRI data include 169 time-point, change it accordingly
unix(['fslcomplex -realcartesian ' directory_m 'cx ' directory_m 'cx_rl ' directory_m 'cx_im 0 168'])

%% Extract motion transformation parameters from magnitude data using FSL MCFlirt
unix(['mcflirt -in ' directory_m 'mg -o ' directory_m 'mg_mcn  -meanvol -spline_final -plots -report -stats -mats'])

%% Use motion transformation parameters to correct real and imaginary data  
directory1_mc = [directory_m 'mg_mcn.mat'];

unix(['applyxfm4D ' directory_m 'cx_rl ' directory_m 'mg_mcn_mean_reg ' directory_m 'cx_rl_mcn ' directory1_mc ' --fourdigit --interp=spline'])
unix(['applyxfm4D ' directory_m 'cx_im ' directory_m 'mg_mcn_mean_reg ' directory_m 'cx_im_mcn ' directory1_mc ' --fourdigit --interp=spline'])

%% Convert motion corrected real and imaginary data into magnitude and phase data  
unix(['fslcomplex -complex ' directory_m 'cx_rl_mcn ' directory_m 'cx_im_mcn ' directory_m 'cx_mc 0 168'])
unix(['fslcomplex -realabs ' directory_m 'cx_mc ' directory_m 'cx_mg_mcn 0 168'])
unix(['fslcomplex -realphase ' directory_m 'cx_mc ' directory_p 'cx_ph_mcn 0 168'])
