function xx = twomey_tikhonov(kk)

xx = ones(1,kk); xm2 = diag(xx,-2);       xm2 = xm2(1:kk,1:kk);
xx = ones(1,kk); xm1 = diag(xx,-1)*(-4);  xm1 = xm1(1:kk,1:kk);
xx = ones(1,kk); x0  = diag(xx,0)*6;      
xx = ones(1,kk); xp1 = diag(xx,1)*(-4);   xp1 = xp1(1:kk,1:kk);
xx = ones(1,kk); xp2 = diag(xx,2);        xp2 = xp2(1:kk,1:kk);

x0(1,1) = 1;   x0(kk,kk) = 1; x0(2,2) = 5; x0(kk-1,kk-1) = 5;
xm1(2,1) = -2; xm1(kk,kk-1) = -2;
xp1(1,2) = -2; xp1(kk-1,kk) = -2;

xx = xm2 + xm1 + x0 + xp1 + xp2;