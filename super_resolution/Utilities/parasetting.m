function par  =   parasetting(par)

par.win           =   4;  % initial patch size: 4*4.
par.rwin          =   4;  % residual patch size: 4*4
par.step          =   3;  % stride
par.rstep         =   2;  % residual stride
% par.every_nblk    =   80;
% par.nblk          =   80;
par.Nsig          =   5;  % variance of the noise image
par.delta         =   0.001;
par.eta           =   2;
par.iter          =   1;
par.f_j           =   10;
par.sigma_win     =   16; % Omage_i: 16

factor            = par.up_scale;
if factor == 2
    par.every_nblk    =   40;
    par.nblk          =   40;
elseif factor == 3
    par.every_nblk    =   80;
    par.nblk          =   80;
else
    par.every_nblk    =   100;
    par.nblk          =   100;   
end

if par.up_scale == 2
   par.num_img       =   2;
else
   par.num_img       =   1;
end

end
