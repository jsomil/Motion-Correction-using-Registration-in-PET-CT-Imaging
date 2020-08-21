%load('ld_CT.mat');
load('nld_CT30.mat');

b=(data1);
A=data;

[x y z] = bitonic(b, 3);

% Show results
for i = 1:3
  
  switch (i)
    case 1
      title = 'Gaussian';
      B = y;
    case 2
      title = 'Median';
      B = z;
    case 3
      title = 'Bitonic';
      B = x;
  end
  
  figure(i);
  imagesc(B);
  colormap('gray');
  axis image
  
end
[peaksnr1, snr] = psnr(y, A);
 mse1=immse(y, A);
 rmse1=sqrt(mse1);
[peaksnr2, snr] = psnr(z, A);
[peaksnr3, snr] = psnr(x, A);
[peaksnr, snr] = psnr(b, A);

[ssimval, ssimmap] = ssim(b,A);

[ssimval1, ssimmap] = ssim(y,A);
[ssimval2, ssimmap] = ssim(z,A);
[ssimval3, ssimmap] = ssim(x,A);

peaksnr
peaksnr1
peaksnr2
peaksnr3

ssimval
ssimval1
ssimval2
ssimval3
rmse1
%%
function out = img_norm(in)
    a = double(in(:));
    a= (a-min(a))./(max(a)-min(a));
    out = reshape(a,size(in));
end




