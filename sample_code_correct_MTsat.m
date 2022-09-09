%% Sample code to correct B1+ inhomogeneity in MTsat maps 
% Please see the README file to make sure you have the necessary MATLAB
% packages to run this code.

% This script is to analyze MTw images obtained at different B1 pulse
% amplitudes applied for the MT saturation pulses. 

% currently set up to be run section by section so that you can view the
% results/ check data quality/ troubleshoot. To run faster, comment out lines used to
% display the images.

% code that can be used to load/export MRI data is here: https://github.com/ulrikls/niak
% image view code is a modified version of https://www.mathworks.com/matlabcentral/fileexchange/47463-imshow3dfull

%% Load the associated matlab files
% Determine where your m-file's folder is.
folder = fileparts(which(mfilename)); 
% Add that folder plus all subfolders to the path.
addpath(genpath(folder));

%% load images
DATADIR = '/folder/where/your/images/are/'; % -> USER DEFINED

% image names
seg_fn = {'high_flip_angle_noMT.mnc.gz' 'low_flip_angle_noMT.mnc.gz' 'MT_weighted_image.mnc.gz'}; % -> USER DEFINED

%load the images your favourite way
% For nifti files use: niak_read_nifti
[hdr, hfa] = niak_read_minc(strcat(DATADIR,seg_fn{1}));
[~, lfa] = niak_read_minc(strcat(DATADIR,seg_fn{2}));
[~, mtw] = niak_read_minc(strcat(DATADIR,seg_fn{3}));

% can check to see if it loaded properly
figure; imshow3Dfull(lfa, [200 600],jet)

%% Load B1 map and set up b1 matrices

% B1 nominal and measured
b1_rms = 2.36; % value in microTesla. Nominal value for the MTsat pulses  % -> USER DEFINED

% load B1 map
[~, b1] = niak_read_minc(strcat(DATADIR,'b1field_filename.mnc'));

% filter the b1 map if you wish. 
b1 = imgaussfilt3(b1,1);

% check the quality of the map, and make sure it is oriented the same as
% your other images (and same matrix size. if not you may need to pad the
% image). 
figure; imshow3Dfull(b1, [0.7 1.2],jet)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% include any preprocessing steps here such as MP-PCA denoising, and
% unringing of GRE images. Not needed for this analysis. Each of the below
% functions can be optimized for your dataset. 
%% denoise (code here https://github.com/sunenj/MP-PCA-Denoising)  % -> USER DEFINED OPTION
% works better with the more images you have. 
img_dn = cat(4,hfa,lfa,mtw);
all_PCAcorr = MPdenoising(img_dn);

%% unring the images ( code here https://github.com/josephdviviano/unring)
% depending on how your images were collected/loaded, the final number in
% the line could be a 1,2 or 3. 
hfa_proc= unring3D(all_PCAcorr(:,:,:,1),3);
lfa_proc= unring3D(all_PCAcorr(:,:,:,2),3);
mtw_proc = unring3D(all_PCAcorr(:,:,:,3) ,3);

% Use the code below to see if you got the right value 
figure; imshow3Dfull(mtw_proc, [150 600])
figure; imshow3Dfull(all_PCAcorr(:,:,:,3), [150 600])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Generate a brain mask to remove background
mask = zeros(size(lfa)); 
mask (lfa >175) = 1;  % check your threshold here, data dependent. You could also load a mask made externally instead. 

% check it 
figure; imshow3Dfullseg(lfa, [150 600], mask)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Begin MTsat calculation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Calculate A0 and R1
low_flip_angle = 5;    % flip angle in degrees % -> USER DEFINED
high_flip_angle = 20;  % flip angle in degrees % -> USER DEFINED
TR = 30;               % repetition time of the GRE kernel in milliseconds % -> USER DEFINED

a1 = low_flip_angle*pi/180 .* b1;
a2 = high_flip_angle*pi/180 .* b1; 

R1 = 0.5 .* (hfa.*a2./ TR - lfa.*a1./TR) ./ (lfa./(a1) - hfa./(a2));
R1 = R1.*mask;
T1 = 1/R1  .* mask;

App = lfa .* hfa .* (TR .* a2./a1 - TR.* a1./a2) ./ (hfa.* TR .*a2 - lfa.* TR .*a1);
App = App .* mask;

%check them
figure; imshow3Dfull(T1, [0 3500],jet)
figure; imshow3Dfull(App , [2500 6000])

% can export your T1 or R1 map here if you wish
% note these results are in milliseconds (T1) or 1/ms (for R1)
hdr.file_name = strcat(DATADIR,'calc_sat_maps/R1.mnc'); niak_write_minc3D(hdr,R1); % -> USER DEFINED

%% Generate MTsat maps for the MTw images. 
% Inital Parameters
readout_flip = 9; % flip angle used in the MTw image, in degrees % -> USER DEFINED
TR = 28; % -> USER DEFINED

% calculate maps as per Helms et al 2008. 
a_MTw_r = readout_flip /180 *pi;
MTsat = (App.* (a_MTw_r*b1)./ mtw_proc - 1) .* (R1) .* TR - ((a_MTw_r*b1).^2)/2;

% check result
figure; imshow3Dfull(MTsat, [0 0.03],jet)

%fix limits for background
MTsat(MTsat<0) = 0;

%% Generate MTsat correction factor maps. 

fitValues = load('fitValuesDirectory/fitValues.mat'); % -> USER DEFINED
fitValues = fitValues.fitValues;

R1_s = R1* 1000; % convert from 1/ms to 1/s

%% Generate MTsat correction factor map. 
CF_MTsat = MTsat_B1corr_factor_map(b1, R1_s, b1_rms,fitValues);

%% Correct the maps
MTsat_b1corr  = MTsat  .* (1+ CF_MTsat)  .* mask;

% display the corrected and uncorrected for comparison
figure; imshow3Dfull(MTsat, [0 0.03],jet); 
figure; imshow3Dfull(MTsat_b1corr, [0 0.03],jet)





