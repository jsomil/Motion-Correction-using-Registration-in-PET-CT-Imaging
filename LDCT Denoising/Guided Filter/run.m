% example: flash/noflash denoising

close all;

%load('ld_CT.mat');
load('nld_CT30.mat');

b=data1;
A=data;% I = double(pet) / 255;
% 
% p = double(ct) / 255;

I = b;
p = b;

r = 8;
eps = 0.02^2;

q = zeros(size(I));
tic;
q(:, :, 1) = guidedfilter(p(:, :, 1), p(:, :, 1), r, eps);
%q(:, :, 2) = guidedfilter(I(:, :, 2), p(:, :, 2), r, eps);
%q(:, :, 3) = guidedfilter(I(:, :, 3), p(:, :, 3), r, eps);
q_sub = q + 22.*(p-q);
toc;
r =8;
s = 4;
% q_sub = zeros(size(I));
% tic;
% q_sub(:, :, 1) = fastguidedfilter(I(:, :, 1), p(:, :, 1), r, eps, s);
% %q_sub(:, :, 2) = fastguidedfilter(I(:, :, 2), p(:, :, 2), r, eps, s);
% %q_sub(:, :, 3) = fastguidedfilter(I(:, :, 3), p(:, :, 3), r, eps, s);
% toc;
%%
figure();
ct=q_sub;
%save('de_pet.mat','pet')
colormap('gray')
imagesc(q);
axis image

[peaksnr1, snr] = psnr(p, A);
[peaksnr2, snr] = psnr(q, A);
mse=immse(q, A);
rmse=sqrt(mse);
[peaksnr3, snr] = psnr(q_sub, A);

[ssimval1, ssimmap] = ssim(p,A);
[ssimval2, ssimmap] = ssim(q,A);
[ssimval3, ssimmap] = ssim(q_sub,A);

peaksnr1
peaksnr2
peaksnr3

ssimval1
ssimval2
ssimval3

rmse

%%
function out = img_norm(in)
    a = double(in(:));
    a= (a-min(a))./(max(a)-min(a));
    out = reshape(a,size(in));
end
