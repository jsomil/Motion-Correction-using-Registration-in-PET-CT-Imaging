theta = linspace(0,360,900);
img = im2double(phantom(512));
output_size=max(size(img));
figure; imshow(img)
c = radon(img, theta);
figure;imshow(c,[]);
a = iradon(c, theta,output_size);
c = 1e12.*imnoise(c./1e12, 'poisson');
figure;imshow(c,[]);
c = iradon(c, theta,output_size);
figure; imshow(c);

data=img_norm(img);
data1=img_norm(c);

save('nld_CT30','data','data1');

RMSE = sqrt(immse(img_norm(a),img_norm(img)));
display(RMSE)
RMSE1 = sqrt(immse(img_norm(c),img_norm(img)));
display(RMSE1)

[peaksnr, snr]=psnr(img_norm(c),img_norm(img));
[ssimval, ssimmap]=ssim(img_norm(c),img_norm(img));

peaksnr
ssimval

%%
function out = img_norm(in)
    a = double(in(:));
    a= (a-min(a))./(max(a)-min(a));
    out = reshape(a,size(in));
end
