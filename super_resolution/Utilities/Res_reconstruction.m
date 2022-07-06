function [im_re] = Res_reconstruction( par )

rwin             = par.rwin;     % patch size
N_patch          = par.N-rwin+1; 
M_patch          = par.M-rwin+1;
Nsig             = par.Nsig;     % noise
im               = zeros(par.N,par.M);
im_wei           = zeros(par.N,par.M);
Xres             = par.Xres;
im_re            = par.res_h; 
y                = par.res_h;
delta            = par.delta;
eta              = par.eta;

for iter = 1 : par.iter
    
    im_re         = im_re + delta*(y - im_re);
    
    if iter == 1
        res_patch   = im_res2patch(im_re,par);
    else
        res_patch   = im_res2patch(im_re,par);
        Xres(:,1:N_patch*M_patch)  = res_patch(:,:);
    end
    
    pos_arr      = match_resblock(res_patch,Xres,par);
    Ys           = zeros(size(res_patch));
    T            = zeros(size(res_patch));
   
    for ii = 1:size(pos_arr,2)
        pos_idx                   =     pos_arr(:,ii);
        pos_idx_self              =     pos_idx(pos_idx <= N_patch * M_patch);
        pos_idx_cooperation       =     pos_idx(pos_idx > N_patch * M_patch);
        RX_self                   =     Xres(:,pos_idx_self);
        RX_cooperation            =     Xres(:,pos_idx_cooperation);
        L                         =     zeros(size(RX_self));
        RX_sc                     =    [RX_self,RX_cooperation];
        if iter == 1
            U{ii} = zeros(size(RX_sc));
        end
        
        [Lii,sigma_L,sigma_RX]    =   solve_NLR(RX_sc+U{ii}/(2*eta),par);
        RRX_sc{ii}                =   RX_sc;
        LL{ii}                    =   Lii;
        L(:,:)                    =   Lii(:,1:size(RX_self,2));
        Ys(:,pos_idx_self) =  Ys(:,pos_idx_self) + L ;
        T(:,pos_idx_self)  =  T(:, pos_idx_self) + 1 ;
    end
    
    for i = 1:rwin
        for j = 1:rwin
            im(i:N_patch+i-1,j:M_patch+j-1)     = im(i:N_patch+i-1,j:M_patch+j-1) + reshape(Ys((i-1)*rwin+j,:)',[N_patch,M_patch]); 
            im_wei(i:N_patch+i-1,j:M_patch+j-1) = im_wei(i:N_patch+i-1,j:M_patch+j-1) + reshape(T((i-1)*rwin+j,:)',[N_patch,M_patch]);
        end
    end

    im_re      = (im_re./(2 *eta* Nsig^2) + im) ./(1/(2 *eta * Nsig^2)+im_wei+eps);

    for ii = 1 : size(pos_arr,2)
            U{ii}    =  U{ii} + eta*(RRX_sc{ii}-LL{ii});
    end
    eta = eta + eta * 0.05;
end

% PINNM
function [Lii,Sigma_L,Sigma_RX]   =   solve_NLR(RX, par)
    [~,~]            =   size(RX);
    [U0,Sigma_0,V0]  =   svd(full(RX),'econ');
    Sigma_RX         =   diag(Sigma_0);
    Sigma_L          =   zeros(size(Sigma_RX));
    f_j              =   par.f_j;
    Sigma_Lj         =   0;
    for iii = 1:size(Sigma_RX)-par.sigma_win+1
        for jjj = iii:iii+par.sigma_win-1
        index_t      =   iii:iii+par.sigma_win-1;
        Sigma_m      =   mean(Sigma_RX(index_t));
        h            =   sum((Sigma_RX(index_t)- Sigma_m).^2)/par.sigma_win;
        h_i          =   sqrt(h);
        Sigma_Lj     =   max((f_j^2 .* Sigma_RX(jjj) + h_i^2 .* Sigma_RX(iii))./(f_j^2+h_i^2), 0) + Sigma_Lj;
        end
        Sigma_L(iii)  = Sigma_Lj ./ par.sigma_win;
    end
     
    rr               =   sum(Sigma_L > 10);
    Lii              =   U0(:,1:rr)*diag(Sigma_L(1:rr))*V0(:,1:rr)';
    
end

end