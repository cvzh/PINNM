function X_sub = im_res2patch(image,par)

win   = par.rwin;   % rwin: the size of residual patch
count = 1;
[height,width] = size(image);
N1 = height-win+1;
M1 = width-win+1;
X_sub = zeros(win^2, N1*M1);

for i = 1:win
    for j = 1:win
        sub_img = image(i:N1+i-1,j:M1+j-1);
        X_sub(count,:)  = sub_img(:)';
        count = count+1;
    end    
end

end
