function [rodgers_rate,errorx,dofs,cdofs,gain,ak,r,se,inv_se,se_errors,ak_water,ak_temp] = rodgers(driver,aux)
%---------------------------------------------------------------------------
% OEM retrieval for RATES so y(x) = sum(rates(i) * jac(:,i)), to compare to yIN
%---------------------------------------------------------------------------
% Notation consistent with 
% Tilman Steck, Methods for determining regularization for atmospheric 
%               retieval problems, Appl Optics, v41, pg 1788 (2002)
%
% some params used in the code
%   jacobian           = k
%   obs error matrix   = se = 2378x2378
%   inverse of this    = diag(1./diag(se)) = inv_se = 2378x2378
%   se_errors          = uncertainties used in building up Se matrix (from obs and forward model)
%                        se_errors.ncerrors = ncerrors;                                       2378x1
%                        se_errors.fmerrors = ones(size(ncerrors)) * driver.oem.sarta_error;  2378x1
%
%   param error matrix = r = 200x200; this starts out as being read in from a mat file 
%                        (typically L0/L1) and then gets manipulated through eg lambdas
%   apriori            = xset = 200x1
%   inital value       = xn   = 200x1
%
% output
%   errorx             = proagated uncertainties in form of a matrix 200x200
%   rodgers_rate       = fitted rates afetr 1 iteration, xnp1 = xn + deltax
%   deg of freedom     = dofs
%   diag(deg freedom)  = cdofs
%   gain matrix        = gain
%   averaging kernel   = ak
%   inds               = actual channels used (maybe slightly different than what user specified, 
%                        if code finds bad rates)
%   r                  = actual relaxation matrix used for parameters (colgas,ST,WV(z),T(z) etc)
%   se                 = actual channel covariance matrix used for spectral obs (eg 1x2378 or 1x8461)
%   inv_se             = actual channel inverse covariance matrix used
%   se_errors          = actual channel uncertainties used
%
%---------------------------------------------------------------------------

% Jacobians
m_ts_jac = aux.m_ts_jac;

% Index of frequencies used
inds     = driver.jacobian.chanset;

[mm,nn] = size(m_ts_jac);
if (length(inds) < nn & driver.oem.dofit)
  fprintf(1,'  length(inds) = %4i \n',length(inds));
  fprintf(1,'  size jac mm,nn = %4i %4i\n',mm,nn);
  disp('Error: More jacobians than channels!')
end

% Apriori state
xb       = aux.xb;

% Covariance (uncertainties/correlations) of measurements
lenr = length(inds);
fme  = ones(1,lenr)*driver.oem.sarta_error;
fme  = diag(fme);          

sizer = size(driver.rateset.rates);
se_errors.fmerrors = ones(sizer) * driver.oem.sarta_error;

%% get 2378x2378 spectral cov matrix
e0 = diag(driver.rateset.unc_rates(inds));
for i=1:length(inds)
   for j=1:length(inds)
      if i~=j
         e0(i,j) = 0.5*sqrt(e0(i,i))*sqrt(e0(j,j));
      end
   end
end
% 
% keyboard

% Error correlation matrix of observations (diagonal)
se = e0 + fme;  
se = se.*se;

% xb is the apriori
[zz1,zz2] = size(xb);
% Linearization point = zero, assuming fits are linear
% note by SSM on 7/4/2013
%   this is a little odd, and makes the code less general purpose!!!!
%   I'd prefer xn = xb!!! of course if xb = 0 this is moot
% xn = zeros(zz1,zz2);  %% orig, before July 2013
xn = xb;              %% after July 2013

% Form k matrix (Jacobians)
k = m_ts_jac(inds,:);
[mm,nn] = size(k);

% Form y - F(xa)
fx = zeros(size(driver.rateset.rates));
for iy = 1 : length(xn)
   fx = fx + (xn(iy)*m_ts_jac(:,iy));
end
deltan = driver.rateset.rates(inds) - fx(inds);

% Do this once to save time, assume diagonal, no need for pinv   ORIG 201
inv_se = pinv(se);          disp(' <<<<<<<< inv_se = pinv(se)');           %%% NEW  post Dec 2012
oo = find(isinf(inv_se) | isnan(inv_se)); inv_se(oo) = 0;

r = inv(driver.oem.cov);
% Use following line for both cov and Tikhonov reg.
% Need to input 1E2 and 1E1 alpha variables via driver.oem
%l = get_l(97,1);s = transpose(l)*l;rc = blkdiag(zeros(6,6),1E2*s,1E1*s);r = r + rc;
% Use following line for only Tikhonov reg.
l = get_l(97,1);
s = transpose(l)*l;
rc = blkdiag(zeros(6,6),driver.oem.alpha_water*s,driver.oem.alpha_temp*s);

switch driver.oem.reg_type
  case 'reg_and_cov'
    r = r + rc;
  case 'reg'
    r = rc;
  case 'cov'
    r = r;
  otherwise
    disp('Incorrect choice driver.oem.reg_type')
end

for ii = 1 : driver.oem.nloop
  % Do the retrieval inversion
  dx1    = r + k' * inv_se * k; 
  dx1    = pinv(dx1);    
  dx2    = k' * inv_se * deltan - r*(xn-xb);
  deltax = dx1*dx2; 

% Update first guess with deltax changes
  rodgers_rate = real(xn + deltax);
  xn = rodgers_rate;

  xsave(ii,:) = rodgers_rate;

  if ii <= driver.oem.nloop
    deltanIN = deltan;
    xn = rodgers_rate;

    % Form the computed rates; also see lines 108-121 of oem_lls.m
    thefitr = zeros(1,length(driver.rateset.rates));
    for ix = 1 : length(xn)
      thefitr = thefitr + xn(ix)*m_ts_jac(:,ix)';
    end

    % Compute chisqr, and new deltan
    deltan = driver.rateset.rates' - thefitr;
    deltan = deltan(:,inds)';
    chisqr(ii) = sum(deltan'.*deltan');
   
    if driver.oem.doplots
       clf
       plot(1:length(deltan),deltanIN,1:length(deltan),deltan,'r'); 
       title(['obs - fit at iteration ' num2str(ii)]); pause(0.1)
   end
  end
end

best = find(chisqr ==  min(chisqr),1);
rodgers_rate = xsave(best,:);

if driver.oem.nloop > 1
  disp('printing out successive chisqr values (upto N-1 th iterate) ...')
  chisqr
end

% Error analysis and diagnostics
errorx = inv(k' * inv_se * k + r);    %% decided pinv is too unstable, Nov 2013, but not much difference
dofsx  = errorx * r; 
dofsx  = eye(size(dofsx)) - dofsx; 
dofs   = trace(dofsx);
cdofs  = diag(dofsx);                 %% so we can do cumulative d.of.f

% Gain is relative weight of first guess and observations
r_water=r(driver.jacobian.water_i,driver.jacobian.water_i); 
r_temp=r(driver.jacobian.temp_i,driver.jacobian.temp_i); 
inv_r = pinv(r);
inv_r_water=pinv(r_water); 
inv_r_temp=pinv(r_temp); 

% inv operator seems OK for this matrix; if problems go back to pinv
k_water=k(:,driver.jacobian.water_i); 
k_temp=k(:,driver.jacobian.temp_i); 
gain  = inv_r *k' * inv(k * inv_r * k' + se);
gain_water=inv_r_water*k_water'*inv(k_water*inv_r_water*k_water'+se); 
gain_temp=inv_r_temp*k_temp'*inv(k_temp*inv_r_temp*k_temp'+se); 

% Compute averaging kernel
ak = gain * k;   
ak_water=gain_water*k_water; 
ak_temp=gain_temp*k_temp; 

