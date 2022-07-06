function par  =   parasetting(par)

par.win           =   4;  % initial patch size: 4*4.
par.rwin          =   4;  % residual patch size: 4*4
par.step          =   3;  % stride
par.rstep         =   2;  % residual stride
% par.every_nblk    =   200;
% par.nblk          =   200;
par.Nsig          =   5;  % variance of the noise image
par.delta         =   0.001;
par.eta           =   2;
par.iter          =   1;
par.f_j           =   10;
par.sigma_win     =   16; % Omage_i: 16

noise_level         =   par.nSig;
if noise_level <= 20
    par.every_nblk    =   25;
    par.nblk          =   25;
elseif noise_level <= 30
    par.every_nblk    =   80;
    par.nblk          =   80;
else
    par.every_nblk    =   160;
    par.nblk          =   160;   
end

if par.nSig  == 15
   par.num_img       =   1;
else
   par.num_img       =   1;
end

end
