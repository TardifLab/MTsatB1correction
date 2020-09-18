%% Sample code to correct B1+ inhomogeneity in MTsat maps from provided sample data

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
DATADIR = '/folder/where/your/images/are/';

% image names
seg_fn = {'hfa.mnc.gz' 'lfa.mnc.gz' 'dual_6p8.mnc.gz' 'dual_7p65.mnc.gz' 'dual_8p5.mnc.gz' 'dual_9p35.mnc.gz' 'pos_6p8.mnc.gz' 'pos_7p65.mnc.gz' 'pos_8p5.mnc.gz' 'irl_2k.mnc.gz'}; 

%load the images your favourite way
[hdr, hfa] = niak_read_minc(strcat(DATADIR,seg_fn{1}));
[~, lfa] = niak_read_minc(strcat(DATADIR,seg_fn{2}));

[~, mtw1_dual] = niak_read_minc(strcat(DATADIR,seg_fn{3}));
[~, mtw2_dual] = niak_read_minc(strcat(DATADIR,seg_fn{4}));
[~, mtw3_dual] = niak_read_minc(strcat(DATADIR,seg_fn{5}));
[~, mtw4_dual] = niak_read_minc(strcat(DATADIR,seg_fn{6}));

[~, mtw1_single] = niak_read_minc(strcat(DATADIR,seg_fn{7}));
[~, mtw2_single] = niak_read_minc(strcat(DATADIR,seg_fn{8}));
[~, mtw3_single] = niak_read_minc(strcat(DATADIR,seg_fn{9}));
[~, mtw4_single] = niak_read_minc(strcat(DATADIR,seg_fn{10}));

[~, mtw_2k] = niak_read_minc(strcat(DATADIR,seg_fn{11}));

% can check to see if it loaded properly, don't worry about orientation
figure; imshow3Dfull(lfa, [200 600],jet)

%% Load B1 map and set up b1 matrices

% B1 nominal and measured
b1_rms = [6.8 7.65 8.5 9.35];  % value in microTesla. Nominal value for the MTsat pulses 

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
%% denoise (code here https://github.com/sunenj/MP-PCA-Denoising) 
img_dn = cat(4,hfa, lfa, mtw1_dual, mtw2_dual, mtw3_dual, mtw4_dual, mtw1_single, mtw2_single, mtw3_single, mtw4_single,mtw_2k);
all_PCAcorr = MPdenoising(img_dn);

%% unring the images ( code here https://github.com/josephdviviano/unring)
% these values are not used as the sample data was already processed
% through these steps. This code is provided simply as reference. 
hfa_proc= unring3D(all_PCAcorr(:,:,:,1),3);
lfa_proc= unring3D(all_PCAcorr(:,:,:,2),3);
mtw1_dual_proc = unring3D(all_PCAcorr(:,:,:,3) ,3);
mtw2_dual_proc = unring3D(all_PCAcorr(:,:,:,4) ,3);
mtw3_dual_proc = unring3D(all_PCAcorr(:,:,:,5) ,3);
mtw4_dual_proc = unring3D(all_PCAcorr(:,:,:,6) ,3);

mtw1_single_proc = unring3D(all_PCAcorr(:,:,:,7) ,3);
mtw2_single_proc = unring3D(all_PCAcorr(:,:,:,8) ,3);
mtw3_single_proc = unring3D(all_PCAcorr(:,:,:,9) ,3);
mtw4_single_proc = unring3D(all_PCAcorr(:,:,:,10) ,3);

mtw_2k_proc = unring3D(all_PCAcorr(:,:,:,11) ,3);
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
low_flip_angle = 5;    % flip angle in degrees
high_flip_angle = 20;  % flip angle in degrees
TR = 30;               % repetition time of the GRE kernel in milliseconds

a1 = low_flip_angle*pi/180 .* b1; % note the inclusion of b1 here.
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
hdr.file_name = strcat(DATADIR,'calc_sat_maps/R1.mnc'); niak_write_minc3D(hdr,R1);

%% Generate MTsat maps for the MTw images. 
% Inital Parameters
readout_flip = 9; % flip angle used in the MTw image, in degrees
TR = 28;
a_MTw_r = readout_flip /180 *pi;


% calculate maps as per Helms et al 2008. Note: b1 is included here for
% flip angle

MTsat_dual1 = (App.* (a_MTw_r*b1)./ mtw1_dual - 1) .* (R1) .* TR - ((a_MTw_r*b1).^2)/2;
MTsat_dual2 = (App.* (a_MTw_r*b1)./ mtw2_dual - 1) .* (R1) .* TR - ((a_MTw_r*b1).^2)/2;
MTsat_dual3 = (App.* (a_MTw_r*b1)./ mtw3_dual - 1) .* (R1) .* TR - ((a_MTw_r*b1).^2)/2;
MTsat_dual4 = (App.* (a_MTw_r*b1)./ mtw4_dual - 1) .* (R1) .* TR - ((a_MTw_r*b1).^2)/2;

MTsat_single1 = (App.* (a_MTw_r*b1)./ mtw1_single - 1) .* (R1) .* TR - ((a_MTw_r*b1).^2)/2;
MTsat_single2 = (App.* (a_MTw_r*b1)./ mtw2_single - 1) .* (R1) .* TR - ((a_MTw_r*b1).^2)/2;
MTsat_single3 = (App.* (a_MTw_r*b1)./ mtw3_single - 1) .* (R1) .* TR - ((a_MTw_r*b1).^2)/2;
MTsat_single4 = (App.* (a_MTw_r*b1)./ mtw4_single - 1) .* (R1) .* TR - ((a_MTw_r*b1).^2)/2;

MTsat_2k = (App.* (a_MTw_r*b1)./ mtw_2k - 1) .* (R1) .* TR - ((a_MTw_r*b1).^2)/2;
% check them, did it work?
%figure; imshow3Dfull(MTsat, [0 0.03],jet)

%fix limits - helps with background noise
MTsat_dual1(MTsat_dual1<0) = 0;
MTsat_dual2(MTsat_dual2<0) = 0;
MTsat_dual3(MTsat_dual3<0) = 0;
MTsat_dual4(MTsat_dual4<0) = 0;

MTsat_single1(MTsat_single1<0) = 0;
MTsat_single2(MTsat_single2<0) = 0;
MTsat_single3(MTsat_single3<0) = 0;
MTsat_single4(MTsat_single4<0) = 0;

MTsat_2k(MTsat_2k<0) = 0;

%% Start with M0B,app fitting  MTsat values 

% Can run the fit with different options, just make sure to load the right
% fit...
dual_s = cat(4, msat_dual_6p8,msat_dual_7p65 ,msat_dual_8p5 , msat_dual_9p35);
pos = cat(4, msat_pos_6p8, msat_pos_7p65, msat_pos_8p5,msat_pos_9p35);

R1_s = R1*1000; % need to convert to 1/s from 1/ms

% load in the fit results from simSeq_M0b_R1obs.m
% makes sure to find the file names/locations where you made the values. 
fitvalsDir = '/FileDirectory/';

fitValues_single = load(strcat(fitvalsDir,'fitValues_7kHz_single.mat'));
fitValues_dual = load(strcat(fitvalsDir,'fitValues_7kHz_dual.mat'));
fitValues2k = load(strcat(fitvalsDir,'fitValues_2kHz.mat'));
fitValues_single = fitValues_single.fitValues;
fitValues_dual = fitValues_dual.fitValues;
fitValues2k = fitValues2k.fitValues;

% initialize matrices
M0b_s = zeros(size(lfa));
M0b_d = zeros(size(lfa));
M0b_2k = zeros(size(lfa));

fit_qual_s = zeros(size(lfa));
fit_qual_d = zeros(size(lfa));
fit_qual_2k = zeros(size(lfa));

comb_res_s = zeros(size(lfa));
comb_res_d = zeros(size(lfa));
comb_res_2k = zeros(size(lfa));

tic % ~ 12 hours to run for 3. 
for i = 1:size(lfa,1)
    for j = 1:size(lfa,2) % for axial slices
        for k = 1:size(lfa,3) % sagital slices
            
            if mask(i,j,k) > 0 %&& dual_s(i,j,k,3) > 0

                [M0b_s(i,j,k), fit_qual_s(i,j,k), comb_res_s(i,j,k)] = CR_fit_M0b_v1( squeeze(b1_comb_scaled(i,j,k,:)), R1_s(i,j,k), squeeze(pos(i,j,k,:)),fitValues_single);
                [M0b_d(i,j,k), fit_qual_d(i,j,k), comb_res_d(i,j,k)] = CR_fit_M0b_v1( squeeze(b1_comb_scaled(i,j,k,:)), R1_s(i,j,k), squeeze(dual_s(i,j,k,:)),fitValues_dual);
                [M0b_2k(i,j,k),fit_qual_2k(i,j,k),comb_res_2k(i,j,k)] = CR_fit_M0b_v1( b1_3p26(i,j,k), R1_s(i,j,k), MTsat_2k(i,j,k),fitValues2k); 
                % If running on one datapoint, your line should look like the M0b_2k line.

            end
        end
    end
    i % output to see it is running. 
end
toc 

figure; imshow3Dfull(M0b_s, [0 0.15],jet)
figure; imshow3Dfull(M0b_d, [0 0.15],jet)
figure; imshow3Dfull(M0b_2k, [0 0.15],jet)


% Save results incase you need to go back, since they take a while to generate!
hdr.file_name = strcat(DATADIR,'calc_sat_maps/M0b_singlefit.mnc'); niak_write_minc3D(hdr,M0b_s);
hdr.file_name = strcat(DATADIR,'calc_sat_maps/M0b_dualfit.mnc'); niak_write_minc3D(hdr,M0b_d);
hdr.file_name = strcat(DATADIR,'calc_sat_maps/M0b_2k.mnc'); niak_write_minc3D(hdr,M0b_2k);

hdr.file_name = strcat(DATADIR,'calc_sat_maps/fit_qual_mask_single.mnc'); niak_write_minc3D(hdr,fit_qual_s);
hdr.file_name = strcat(DATADIR,'calc_sat_maps/fit_qual_mask_dual.mnc'); niak_write_minc3D(hdr,fit_qual_d);
hdr.file_name = strcat(DATADIR,'calc_sat_maps/fit_qual_mask_2k.mnc'); niak_write_minc3D(hdr,fit_qual_2k);

hdr.file_name = strcat(DATADIR,'calc_sat_maps/comb_residuals_single.mnc'); niak_write_minc3D(hdr,comb_res_s);
hdr.file_name = strcat(DATADIR,'calc_sat_maps/comb_residuals_dual.mnc'); niak_write_minc3D(hdr,comb_res_d);
hdr.file_name = strcat(DATADIR,'calc_sat_maps/comb_residuals_2k.mnc'); niak_write_minc3D(hdr,comb_res_2k);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Now make some plots of my fitted maps to see how they correlate with R1 values

R1_p = R1(mask>0);
M0b_d_p = M0b_d(mask>0);
M0b_s_p = M0b_s(mask>0);
M0b_2k_p = M0b_2k(mask>0);

plot_con = cat(2, R1_p,M0b_d_p,M0b_s_p,M0b_2k_p); 
contrast_fit = zeros(3,2);

ft = fittype('poly1');

% Generate quality Plots

%M0b from dual
    tmp = plot_con(:,2);
    tmp_r1 = plot_con(:,1)*1000; % convert from ms to sec
    tmp_r1(tmp==0) = []; % to clean up the plots
    tmp(tmp==0) = [];
    M0b_d_fit = fit(tmp_r1,tmp,ft);
    [R,P]= corrcoef([tmp, tmp_r1],'Alpha',0.05,'Rows','complete') ;
    contrast_fit(1,1) = R(2,1);
    contrast_fit(1,2) = P(2,1);
    fitvals_dual = coeffvalues(M0b_d_fit);
        
    figure;
    heatscatter(tmp_r1,tmp, 1,'.')
    xlim([0.35 1.5])
    ylim([0 0.16])
    hold on
    plot(M0b_d_fit,'fit',0.95);
    caption = sprintf('M_{0,app}^B = %.2g * R_1 %.2g', fitvals_dual(1), fitvals_dual(2));
    text(0.45, 0.149, caption, 'FontSize', 16);    
    caption2 = sprintf('r = %.2f', contrast_fit(1,1));
    text(0.78, 0.135, caption2, 'FontSize', 16);
    ax = gca;
    ax.FontSize = 20; 
    xlabel('R_1 (1/s)', 'FontSize', 20, 'FontWeight', 'bold')
    ylabel('M_{0,app}^B', 'FontSize', 20, 'FontWeight', 'bold')
    %colorbar('off')
    legend('hide')
    
    
    
%M0b from single sat
    tmp = plot_con(:,3);
    tmp_r1 = plot_con(:,1) *1000; % convert from ms to sec
    tmp_r1(tmp==0) = []; % to clean up the plots
    tmp(tmp==0) = [];
    M0bs_fit = fit(tmp_r1,tmp,ft)
    [R,P]= corrcoef([tmp, tmp_r1],'Alpha',0.05,'Rows','complete') ;
    contrast_fit(2,1) = R(2,1);
    contrast_fit(2,2) = P(2,1);
    fitvals_single = coeffvalues(M0bs_fit);
        
    figure;
    heatscatter(tmp_r1,tmp, 1,'.')
    ylim([0 0.2])
    xlim([0.35 1.5])
    hold on
    plot(M0bs_fit,'fit',0.95);
    caption = sprintf('M_{0,app}^B = %.2g * R_1 %.2g', fitvals_single(1), fitvals_single(2));
    text(0.4, 0.185, caption, 'FontSize', 16);   
    caption2 = sprintf('r = %.2f', contrast_fit(2,1));
    text(0.78, 0.165, caption2, 'FontSize', 16);
    ax = gca;
    ax.FontSize = 20; 
    xlabel('R_1 (1/s)', 'FontSize', 20, 'FontWeight', 'bold')
    ylabel('M_{0,app}^B', 'FontSize', 20, 'FontWeight', 'bold')
     %   colorbar('off')
    legend('hide')
    
 %M0b from 2k sat
    tmp = plot_con(:,4);
    tmp_r1 = plot_con(:,1) *1000; % convert from ms to sec
    tmp_r1(tmp==0) = []; % to clean up the plots
    tmp(tmp==0) = [];
    M0bs_fit = fit(tmp_r1,tmp,ft)
    [R,P]= corrcoef([tmp, tmp_r1],'Alpha',0.05,'Rows','complete') ;
    contrast_fit(2,1) = R(2,1);
    contrast_fit(2,2) = P(2,1);
    fitvals2k = coeffvalues(M0bs_fit);
    
    figure;
    heatscatter(tmp_r1,tmp, 1,'.')
    ylim([0 0.13])
    xlim([0.35 1.5])
    hold on
    plot(M0bs_fit,'fit',0.95);
    caption = sprintf('M_{0,app}^B = %.2g * R_1 %.2g', fitvals2k(1), fitvals2k(2));
    text(0.4, 0.12, caption, 'FontSize', 16);   
    caption2 = sprintf('r = %.2f', contrast_fit(2,1));
    text(0.78, 0.105, caption2, 'FontSize', 16);
    ax = gca;
    ax.FontSize = 20; 
    xlabel('R_1 (1/s)', 'FontSize', 20, 'FontWeight', 'bold')
    ylabel('M_{0,app}^B', 'FontSize', 20, 'FontWeight', 'bold')
     %   colorbar('off')
    legend('hide')
    
    
%% Now add these regression equations to the fitValues structure and save. 

fitValues_dual.Est_M0b_from_R1 = strcat( num2str(fitvals_dual(1)),' *Raobs + ',num2str(fitvals_dual(2)));
fitValues_single.Est_M0b_from_R1 = strcat( num2str(fitvals_single(1)),' *Raobs + ',num2str(fitvals_single(2))); 
fitValues2k.Est_M0b_from_R1 = strcat( num2str(fitvals2k(1)),' *Raobs + ',num2str(fitvals2k(2)));

save(strcat(fitvalsDir,'fitValues_7kHz_dual.mat'),'fitValues_dual')
save(strcat(fitvalsDir,'fitValues_7kHz_single.mat'),'fitValues_single')
save(strcat(fitvalsDir,'fitValues_2kHz.mat'),'fitValues2k')


















