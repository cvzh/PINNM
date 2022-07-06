
clear; clc; close all;
addpath('Data','Utilities');

nSig = 25; % 15/25/50
par.nSig = nSig;
folder = 'Data/CBSD68'; 
filepaths = dir(fullfile(folder, '*.jpg'));
par.img   = length(filepaths);
par = parasetting(par); 

for n = 1 : par.img
    X       = [];
    Xres    = [];
    im_hr_y = [];
    for m = 1 : par.num_img
        if m == 1
           O_img = imread(fullfile(folder, filepaths(n).name));
        else
           O_img = im_hr_y;
        end
    [re_height,re_width,re_dim] = size(O_img);
    im_re = zeros(re_height, re_width, re_dim);
    if size(O_img,3) > 1
    im_ycbcr   = double(rgb2ycbcr(O_img));
    im_ycbcr_y = im_ycbcr(:,:,1);
    else
    im_ycbcr_y = double(O_img);
    end
    [height,width] = size(im_ycbcr_y);
    
    if m == 1
        par.image_o          = im_ycbcr_y;
        par.image_o_shrink   =  par.image_o + par.nSig*randn(size(par.image_o));
    else
        par.image_o_shrink   =  par.image_o + par.nSig*randn(size(par.image_o));
    end
    
    if m == 1
        par.y = adpmedft(par.image_o_shrink, 19);
        par.N = height;
        par.M = width;
    else
        par.y = im_ycbcr_y;
        par.N = height;
        par.M = width;  
    end
    
    X_sub = im2patch(im_ycbcr_y, par);
    X     = [X, X_sub];
    par.X = X;

    im_re(:,:,1) = SR_reconstruction(par);

   %----------------------Residual_SR----------------------
    im_re_y      = im_re(:,:,1);
    im_re_y_lr   = adpmedft(im_re_y, 19);
    imo = imread(fullfile(folder, filepaths(n).name));
    if size(imo,3) > 1
      imo_ycbcr   = double(rgb2ycbcr(imo));
      imo_y       = imo_ycbcr(:,:,1);
    else
      imo_y       = double(imo);
    end
   
    res       = par.image_o_shrink - im_re_y_lr;
    par.res   = res;
    if m == 1
       res_h     = adpmedft(res, 19);
    else
       im_lr     = adpmedft(res, 19);
       res_h     = im_ycbcr_y - im_lr;
    end
   
    par.res_h = res_h;
    Xres_sub  = im_res2patch(res_h, par);
    Xres      = [Xres, Xres_sub];
    par.Xres  = Xres;

    im_res  = Res_reconstruction(par);
   
    im_hr_y = im_re_y + im_res;
   
    psnr(n) = compute_psnr(imo_y, im_hr_y);
    ssim(n) = compute_ssim(imo_y, im_hr_y);
    disp(filepaths(n).name);
    fprintf('HR_PSNR is : %.4f  \t', psnr(n));
    fprintf('HR_SSIM is : %.4f  \t', ssim(n));
    fprintf('HR_FSIM is : %.4f  \n', FeatureSIM(imo_y,im_hr_y));
    impath = filepaths(n).name;
    imwrite(im_hr_y/255, impath);  
    end
end

fprintf('Ave_PSNR is : %.4f  \t', sum(psnr)/length(psnr));
fprintf('Ave_SSIM is : %.4f  \n', sum(ssim)/length(ssim));
   


