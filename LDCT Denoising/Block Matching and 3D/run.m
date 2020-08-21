% Grayscale BM3D denoising demo file, based on Y. M??kinen, L. Azzari, A. Foi, 2019.
% Exact Transform-Domain Noise Variance for Collaborative Filtering of Stationary Correlated Noise.
% In IEEE International Conference on Image Processing (ICIP), pp. 185-189

% ---
% The location of the BM3D files -- this folder only contains demo data
addpath('bm3d');

% Experiment specifications   
%imagename = 'cameraman256.png';

% Load noise-free image
y = double(dicomread('47530905'));
% load('nld_CT30.mat');
% y=data;
% Possible noise types to be generated 'g0', 'g1', 'g2', 'g3', 'g4', 'g1w',
% 'g2w', 'g3w', 'g4w'.
noise_type =  'g1';
%g1,g4w
noise_var = 0.02; % Noise variance
seed = 0; % seed for pseudorandom noise realization
% Generate noise with given PSD
[noise, PSD, kernel] = getExperimentNoise(noise_type, noise_var, seed, size(y));
% N.B.: For the sake of simulating a more realistic acquisition scenario,
% the generated noise is *not* circulant. Therefore there is a slight
% discrepancy between PSD and the actual PSD computed from infinitely many
% realizations of this noise with different seeds.

% Generate noisy image corrupted by additive spatially correlated noise
% with noise power spectrum PSD

z =  data1;

% Call BM3D With the default settings.
y_est = BM3D(z, PSD);

% To include refiltering:
%y_est = BM3D(z, PSD, 'refilter');

% For other settings, use BM3DProfile.
% profile = BM3DProfile(); % equivalent to profile = BM3DProfile('np');
% profile.gamma = 6;  % redefine value of gamma parameter
% y_est = BM3D(z, PSD, profile);

% Note: For white noise, you may instead of the PSD 
% also pass a standard deviation 
% y_est = BM3D(z, sqrt(noise_var));


        figure();
        colormap('gray');
        imagesc(y); %title("High Dose CT");
        axis image;




        figure();
        colormap('gray');
        imagesc(z); %title("Noisy Image");
        axis image;
        
        figure();
        colormap('gray');
        imagesc(y_est); %title(" Denoised BM3D");
        axis image;



        [peaksnr1, snr] = psnr(y_est, y);
        [ssimval1, ssimmap] = ssim(y_est,y);
        mse1=immse(y_est, y);
        rmse1=sqrt(mse1);
        [peaksnr, snr] = psnr(z, y);
        [ssimval, ssimmap] = ssim(z,y);
        mse=immse(z, y);
        rmse=sqrt(mse);
        peaksnr
        ssimval
        rmse
        peaksnr1
        ssimval1
        rmse1
%%
function out = img_norm(in)
    a = double(in(:));
    a= (a-min(a))./(max(a)-min(a));
    out = reshape(a,size(in));
end



