function operator = derivative_operators(kk,order)

%% Tilman Steck, Methods for determining regularization for atmospheric 
%%               retieval problems, Appl Optics, v41, pg 1788 (2002)

if order == 0
  kkx = kk;
  xx = eye(kkx);
  operator = xx;
elseif order == 1
  kkx = kk+1;
  xx = ones(1,kkx);   xm1 = diag(xx,0) * -1;
  xx = ones(1,kkx-1); xp1 = diag(xx,+1)* +1;
  xfinal = xm1 + xp1;
  xfinal = xfinal(1:kk-1,1:kkx-1);
  xfinal = xfinal' * xfinal;
  operator = xfinal;
elseif order == 2
  kkx = kk+2;
  xx = ones(1,kkx);   x0  = diag(xx,0) * +1;
  xx = ones(1,kkx-1); xp1 = diag(xx,+1)* -2;
  xx = ones(1,kkx-2); xp2 = diag(xx,+2)* +1;
  xfinal = x0 + xp1 + xp2;
  xfinal = xfinal(1:kk-2,1:kkx-2);
  xfinal = xfinal' * xfinal;
  operator = xfinal;
else
  error('can only make first/second order matrices')
end

