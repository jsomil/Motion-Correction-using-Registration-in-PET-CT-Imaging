load('nld_CT30.mat');
img_n=(data1);
%img=double(dicomread('47530905'));
img=data;
k=wiener2(img_n,[3 3]);
%k=medfilt2(img_n,[3 3]);
%k=fft2(img_n);

        figure();
        colormap('gray');
        imagesc(img); %title("high dose");
        axis image;




        figure();
        colormap('gray');
        imagesc(img_n); %title("Noisy");
        axis image;



        [peaksnr1, snr] = psnr(img_n, img);
        [ssimval1, ssimmap] = ssim(img_n,img);
        
        [peaksnr2, snr] = psnr(k, img);
        [ssimval2, ssimmap] = ssim(k,img);
        mse=immse(k, img);
        rmse=sqrt(mse);

        peaksnr1
        peaksnr2
        rmse

        ssimval1
        ssimval2
        
        figure();
        colormap('gray');
        imagesc(k); %title("denoised");
        axis image;
        
  %%
function out = img_norm(in)
    a = double(in(:));
    a= (a-min(a))./(max(a)-min(a));
    out = reshape(a,size(in));
end



