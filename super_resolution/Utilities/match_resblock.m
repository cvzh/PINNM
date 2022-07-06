function pos_arr0 = match_resblock(image_res_patch,Xres,par)

par         = parasetting(par);
m           = par.every_nblk; 
all_nblk    = par.nblk;
rwin        = par.rwin;
step        = par.rstep;
N           = par.N-rwin+1;
M           = par.M-rwin+1;
r           = 1:step:N;
r           = [r r(end)+1:N];
c           = 1:step:M;
c           = [c c(end)+1:M];
S           = 15; 

Xres        =   Xres';
mX          =   sum(Xres.^2, 2)/2;
x           =   image_res_patch';
mx          =   sum(x.^2,2)/2; 

% Index image
I         =   (1:N*M*par.num_img);
I         =   reshape(reshape(I, N*M, par.num_img), N, M, par.num_img);
N1        =   length(r);
M1        =   length(c);

pos_arr   =   []; % Position
dis_arr   =   []; % Distance

for n = 1 
   
   pos_subarr   =   zeros(m, N1*M1);
   dis_subarr   =   zeros(m, N1*M1);

   t = 1;
   for i = 1 : N1
       for j = 1 : M1
        
        off    =   (c(j)-1)*N + r(i) + (n-1) * N * M; 
        rmin    =   max( r(i) - S, 1 );
        rmax    =   min( r(i) + S, N );
        cmin    =   max( c(j) - S, 1 );
        cmax    =   min( c(j) + S, M );
        
        I_sub   =   I(rmin:rmax, cmin:cmax, n);
        idx     =   reshape(I_sub, [], 1);

        if M*N*(n-1) < off && off <= M*N*n
            ii      =   idx==off;
            idx(ii) =   [];
            dis     =   (mX(idx, :) + mx(off - (n-1) * N * M, :) - Xres(idx, :)*x(off - (n-1) * N * M, :)'); 
            [val,ind]               =   sort(dis);
            pos_subarr(2:m, t)      =   idx(ind(1:m-1));
            pos_subarr(1, t)        =   off;
            dis_subarr(2:m, t)      =   val(1:m-1);
            dis_subarr(1, t)        =   0;
        else
            dis     =   (mX(idx, :) + mx(off - (n-1) * N * M, :) - Xres(idx, :)*x(off - (n-1) * N * M, :)'); 
            [val,ind]               =   sort(dis);      
            pos_subarr(1:m, t)      =   idx(ind(1:m));
            dis_subarr(1:m, t)      =   val(1:m);
        end
        t = t + 1;
       end
   end
   
    pos_arr = [pos_arr;pos_subarr];
    dis_arr = [dis_arr;dis_subarr];

end

pos_arr0 = zeros(all_nblk, N1*M1);

for i = 1:size(pos_arr,2)
   [~,ind]            =   sort(dis_arr(:,i));     
   pos_arr0(:,i)      =   pos_arr(ind(1:all_nblk),i);
end

end

