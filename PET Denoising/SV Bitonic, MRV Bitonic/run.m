
% load('pet_reg.mat');
% pet=Ireg;

%load('./PET_CT_Data/ctdata.mat');
%image=ctdata(:,:,293);
%ct = imcrop(image,[193 152 130 160]);
%ct = imcrop(image,[183 142 150 180]);
%ct=imresize(ct,size(pet));

% read in image
%A = imread('varying_noise.png','PNG');
%A = imread('svbitonic2.png','PNG');

load('nld_CT30.mat');

%load('maskedct.mat');

A=data1;
%b=double(dicomread('47530905'));
b=data;
% Run filters and show results
for i = 1:4
  
  switch (i)
    case 1
      title = 'Original image';
      B = A;
      
    case 2
      title = 'Bitonic, f=6';
      fprintf(1,'Calculating bitonic ... ');
      tic;
      B = bitonic2(A, 6);
      x=B;
      t = toc;
      fprintf(1,'in %.3f secs.\n', t);
      
    case 3
      title = 'Structurally varying bitonic, f=10';
      fprintf(1,'Calculating structurally varying bitonic ... ');
      tic;
      B = svbitonic2(A, 10);
      y=B;
      t = toc;
      fprintf(1,'in %.3f secs.\n', t);

    case 4
      title = 'Multi-resolution varying bitonic, f=9';
      fprintf(1,'Calculating multi-resolution varying bitonic ... ');
      tic;
      B = mvbitonic2(A, 9);
      z=B;
      t = toc;
      fprintf(1,'in %.3f secs.\n', t);

  end
  
%   [peaksnr, snr] = psnr(B, b);
%   [ssimval, ssimmap] = ssim(B,b);
%   mse=immse(B,b);
%   rmse=sqrt(mse);
% 
%   peaksnr
%   ssimval
%   rmse
  figure(i);
  imagesc(B);
  colormap(gray);
  axis image;

end

figure();
imagesc(x-A);
colormap(gray);
axis image;

figure();
imagesc(y-A);
colormap(gray);
axis image;


figure();
imagesc(z-A);
colormap(gray);
axis image;

save ct_bitonic.mat x y z
norm(x-A)
norm(y-A)
norm(z-A)



        [peaksnr1, snr] = psnr(y, b);
        [ssimval1, ssimmap] = ssim(y, b);
        mse1=immse(y, b);
        rmse1=sqrt(mse1);
        [peaksnr2, snr] = psnr(z, b);
        [ssimval2, ssimmap] = ssim(z, b);
        mse2=immse(z, b);
        rmse2=sqrt(mse2);
        [peaksnr, snr] = psnr(x, b);
        [ssimval, ssimmap] = ssim(x, b);
        mse=immse(x, b);
        rmse=sqrt(mse);
        
        peaksnr
        peaksnr1
        peaksnr2

        ssimval
        ssimval1
        ssimval2
        
        rmse
        rmse1
        rmse2
%%
function out = img_norm(in)
    a = double(in(:));
    %a= ((a-min(a))./(max(a)-min(a))).*255;
    a= (a-min(a))./(max(a)-min(a));
    out = reshape(a,size(in));
end

