%% Physiological noise correction based on anatomical CompCor
% V Malekian, FIL Physics, March 2023
% "fmri_compcor" and "fmri_cleaning" functions are available at https://github.com/dmascali/fmri_denoising

addpath('/export/home/vmalekian/NIfTI_20140122/')
addpath('/export/home/vmalekian/Matlab_lib/spm12/')
addpath('/export/home/vmalekian/physio_correction/fmri_denoising-master/')

%% Steps 1 & 2
%% Data preparation and pre-processing (FEAT preprocessing using design.fsf) 

%%Motion correction, HP filtering and generate mean fMRI voulme
dir_feat = '/export/home/vmalekian/physio_correction/example/';
unix(['feat ' (dir_feat) 'design.fsf'])

%% Steps 3 & 4
%% Registeration
dir_fmri = '/export/home/vmalekian/physio_correction/example/fMRI.feat/';

unix(['flirt -in ' dir_fmri 'example_func -ref ' dir_feat 'MT_3DEPI -out ' dir_fmri 'func2MTEPI -omat ' dir_fmri 'func2MTEPI.mat -dof 6 -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -interp spline']) 
unix(['convert_xfm -inverse -omat ' dir_fmri 'MTEPI2func.mat ' dir_fmri 'func2MTEPI.mat'])
unix(['bet2 ' dir_fmri 'mean_func ' dir_fmri 'mean_func_brain -f 0.25 -n -m'])

%% Step 5 
%% MT-3DEPI segmentation

matlabbatch{1}.spm.spatial.preproc.channel.vols = {
    [dir_feat 'MT_3DEPI.nii,1']
    }

matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.01;
matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 30;
matlabbatch{1}.spm.spatial.preproc.channel.write = [1 1];
matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = {'/export/home/vmalekian/Matlab_lib/spm12/tpm/TPM.nii,1'};
matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = {'/export/home/vmalekian/Matlab_lib/spm12/tpm/TPM.nii,2'};
matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = {'/export/home/vmalekian/Matlab_lib/spm12/tpm/TPM.nii,3'};
matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = {'/export/home/vmalekian/Matlab_lib/spm12/tpm/TPM.nii,4'};
matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = {'/export/home/vmalekian/Matlab_lib/spm12/tpm/TPM.nii,5'};
matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = {'/export/home/vmalekian/Matlab_lib/spm12/tpm/TPM.nii,6'};
matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 2;
matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.0005 0.25 0.025 0.1];
matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
matlabbatch{1}.spm.spatial.preproc.warp.samp = 1;
matlabbatch{1}.spm.spatial.preproc.warp.write = [0 0];
matlabbatch{1}.spm.spatial.preproc.warp.vox = NaN;
matlabbatch{1}.spm.spatial.preproc.warp.bb = [NaN NaN NaN
    NaN NaN NaN];

spm_jobman('run', matlabbatch);

unix(['fslmaths ' dir_feat 'c2MT_3DEPI -thr 0.99 -bin ' dir_feat 'c2MT_3DEPI_mask']) %% WM
unix(['fslmaths ' dir_feat 'c3MT_3DEPI -thr 0.99 -bin ' dir_feat 'c3MT_3DEPI_mask']) %% CSF

%% Step 6
%% WM anc CSF mask creation and  transformation into the fMRI space

unix(['applywarp --ref=' dir_fmri 'example_func --in=' dir_feat 'c2MT_3DEPI_mask --out=' dir_feat 'c2MT_3DEPI_mask_2fMREPI --premat=' dir_fmri 'MTEPI2func.mat --interp=nn'])
unix(['applywarp --ref=' dir_fmri 'example_func --in=' dir_feat 'c3MT_3DEPI_mask --out=' dir_feat 'c3MT_3DEPI_mask_2fMREPI --premat=' dir_fmri 'MTEPI2func.mat --interp=nn'])

unix(['fslmaths ' dir_feat 'c2MT_3DEPI_mask_2fMREPI -mul ' dir_fmri 'mean_func_brain_mask ' dir_feat 'c2MT_3DEPI_mask_2fMREPI_BET'])
unix(['fslmaths ' dir_feat 'c3MT_3DEPI_mask_2fMREPI -mul ' dir_fmri 'mean_func_brain_mask ' dir_feat 'c3MT_3DEPI_mask_2fMREPI_BET'])

string = 'c2MT_3DEPI_mask_2fMREPI_BET.nii.gz';
a1m=dir([dir_feat,string]);
nf=load_untouch_nii([dir_feat,a1m(1).name]);
roi_white= (nf.img);

string = 'c3MT_3DEPI_mask_2fMREPI_BET.nii.gz';
a1m=dir([dir_feat,string]);
nf=load_untouch_nii([dir_feat,a1m(1).name]);
roi_csf= (nf.img);

%%remove very top and bottom slices
roi_white(:,:,[1:6,end-6:end])= 0;
string_compcor = [dir_feat 'c2MT_3DEPI_mask_2fMREPI_BET_roi.nii.gz'];
save_nii(make_nii((roi_white)),string_compcor)
unix(['fslcpgeom ' dir_fmri 'mean_func_brain_mask.nii.gz ' string_compcor])

roi_csf(:,:,[1:6,end-6:end])= 0;
string_compcor = [dir_feat 'c3MT_3DEPI_mask_2fMREPI_BET_roi.nii.gz'];
save_nii(make_nii((roi_csf)),string_compcor)
unix(['fslcpgeom ' dir_fmri 'mean_func_brain_mask.nii.gz ' string_compcor])

%% Steps 7 & 8 
%% Extracting regressors 

string = 'c2MT_3DEPI_mask_2fMREPI_BET_roi.nii.gz';
a1m=dir([dir_feat,string]);
nf=load_untouch_nii([dir_feat,a1m(1).name]);
roi_W= (nf.img);

string = 'c2MT_3DEPI_mask_2fMREPI_BET_roi.nii.gz';
a1m=dir([dir_feat,string]);
nf=load_untouch_nii([dir_feat,a1m(1).name]);
roi_C= (nf.img);

string = 'filtered_func_data.nii.gz';
a1m=dir([dir_fmri,string]);
nf=load_untouch_nii([dir_fmri,a1m(1).name]);
data= (nf.img);

C_comp =5;
rois = {roi_C,roi_W};
%%Caluclate number of WM regressors
W_comp = fmri_compcor_wm(data,rois,[C_comp C_comp]);
%%Caluclate number of WM regressors
X= fmri_compcor(data,rois,[C_comp W_comp]);% Alternitvely can use X= fmri_compcor(data,rois,[0.25 0.25])
%%ploting 10 CSF and WM regressors
figure,subplot(5,2,1),plot(X(:,1)),subplot(5,2,2),plot(X(:,2)),subplot(5,2,3),plot(X(:,3)),subplot(5,2,4),plot(X(:,4)),subplot(5,2,5),plot(X(:,5))
subplot(5,2,6),plot(X(:,6)),subplot(5,2,7),plot(X(:,7)),subplot(5,2,8),plot(X(:,8)),subplot(5,2,9),plot(X(:,9)),subplot(5,2,10),plot(X(:,10))

%% Step 9
%% Correlation 
Corr_reg = (tril(corr(X)) - eye(size(X,2)))>0.5;
[ic,jc]= find(Corr_reg==1);
X_n= X;
X_n(:,ic)=[];
save('phyiso.txt','X_n')

%% Step 10
%% fMRI data cleaning 

[res,model_info] = fmri_cleaning(data,-1,X_n);

string_compcor = [dir_fmri 'filtered_func_data_cleaned.nii.gz'];
save_nii(make_nii((res)),string_compcor)
unix(['fslcpgeom ' dir_fmri 'filtered_func_data.nii.gz ' string_compcor])
%% Step 11
%% Effective tSNR evaluation

unix(['fslmaths ' dir_fmri 'filtered_func_data -Tmean ' dir_fmri 'filtered_func_data_mean'])
unix(['fslmaths ' dir_fmri 'filtered_func_data -Tstd ' dir_fmri 'filtered_func_data_std'])
unix(['fslmaths ' dir_fmri 'filtered_func_data_mean -div ' dir_fmri 'filtered_func_data_std ' dir_fmri 'filtered_func_data_snr'])


unix(['fslmaths ' dir_fmri 'filtered_func_data_cleaned -Tmean ' dir_fmri 'filtered_func_data_cleaned_mean'])
unix(['fslmaths ' dir_fmri 'filtered_func_data_cleaned -Tstd ' dir_fmri 'filtered_func_data_cleaned_std'])
unix(['fslmaths ' dir_fmri 'filtered_func_data_cleaned_mean -div ' dir_fmri 'filtered_func_data_cleaned_std ' dir_fmri 'filtered_func_data_cleaned_std_snr'])
eff = sqrt(model_info.dof/model_info.N);
unix(['fslmaths ' dir_fmri 'filtered_func_data_cleaned_std_snr -mul ' sprintf('%02d',eff) '  ' dir_fmri 'filtered_func_data_cleaned_std_fsnr'])

