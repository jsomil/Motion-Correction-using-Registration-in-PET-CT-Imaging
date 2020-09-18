%Kolbeinn Karlsson 06/11/12
%Advanced Multimedia Processing (AMP) Lab, Cornell University
clear; close all; clc;

%Time the algorithm
ticID = tic;

%%
% PCT Parameters
PRE = 10;       %Pre-enhancement cutoff (first frame included in calc's)
POST = 145; %89    %Post-enhancement cutoff (last frame included in calc's)
kappa = 0.73;   %Hematocrit correction factor
loth = 0;       %Lower segmentation threshold
hith = 120;     %Upper segmentation threshold
rho = 1.05;     %Average brain tissue density
fsize = 3;      %Size of Gaussian filter (spatial filtering)
sigma = 0.9;    %Standard deviation of Gaussian filter (spatial filtering)
k = 1;          %Contrast conversion factor
ftsize = 5;    %Size of temporal gaussian filter.
dt = 1; %0.5;       %The time interval between CT samples

% SVD parameters
lambda = 0.3;   %Truncation parameter
m = 2;          %Extend the data matrix m time for block circulant

% Min and max values
min_cbv = 0;
max_cbv = Inf;
min_cbf = 0;
max_cbf = Inf;
min_mtt = 0;
max_mtt = Inf;
min_ttp = 0;
max_ttp = Inf;
min_preprocessed = 0;
max_preprocessed = Inf;
min_bbbp = 0;
max_bbbp = Inf;
min_aif = 1;
max_aif = Inf;
max_delay = 15; % average circulation time 6-11 sec

first_bbbp = 73;
last_bbbp = 89;

I0 = 190; % High-dose tube current level
I = 50; % Simulated low-dose tube current level

% y1 = 200; y2 = 202; x1 = 200; x2 = 202;

% Get the data
% load 'Data/IRB_001_S3.mat'
%load '../DICOM/IRB/IRB_017.mat'
%This file contains the following variables:
% V  - containing the image data [T x Y x X]
% aif_x - the x coordinate of the AIF
% aif_y - the y coordinate of the AIF
% vof_x - the x coordinate of the VOF
% vof_y - the y coordinate of the VOF

data =double(dicomread('47530905'));
%aif_x = AIFx; aif_y = AIFy; vof_x = VOFx; vof_y = VOFy;

figure();
colormap('gray');
imagesc(data); title("CT Image");
axis image;


%Get length
[len,h,w] = size(data);

%Simulate low-dose
sigma = pct_mA2sigma(I,I0);
%sigma     - Standard deviation of the noise
%spectral noise
% figure();
% colormap('gray');
% imagesc(sigma); title("Sigma Image CT");
% axis image;

load('acf.mat');
data1 = pct_noise(data,acf,sigma);
% add acf
%respository spd

figure();
colormap('gray');
imagesc(data1); title("Spectral Noise added CT Image");
axis image;

save ld_CT.mat data1

[peaksnr1, snr] = psnr(data1, data);

[ssimval1, ssimmap] = ssim(data1,data);

peaksnr1


ssimval1




figure();
colormap('gray');
imagesc(data1-double(data)); title("added Noise in CT");
axis image;

b=imread('ct.png');
c=imread('noisy_ct.png');

[peaksnr2, snr] = psnr(c, b);

[ssimval2, ssimmap] = ssim(c,b);

peaksnr2
ssimval2
    

%Get time performance
toc(ticID);
clear ticID

