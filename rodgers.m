function [rodgers_rate,errorx,dofs,cdofs,gain,ak,r,se,inv_se,se_errors,ak_water,ak_temp,ak_ozone] = rodgers(driver,aux)

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

invtype = 0;  %% inv
invtype = 1;  %% pinv     BEST *****
invtype = 2;  %% S. Rump      invillco   addpath /home/sergio/MATLABCODE/IntLab
invtype = 3;  %% T. A. Davis  factorize  addpath /home/sergio/MATLABCODE/FactorizeMatrix/Factorize

if ~isfield(aux,'invtype')
  aux.invtype = 1;   %% default is to use pinv
end
invtype = aux.invtype;
if invtype < 0 | invtype > 3
  error('need invtype between 0 and 3')
end
fprintf(1,'inverse of matrices using method (0) inv (1) PINV (default) (2) invillco (3) factorize : %2i \n',invtype);

if invtype == 2
  addpath /home/sergio/MATLABCODE/IntLab
elseif invtype == 3
  addpath /home/sergio/MATLABCODE/FactorizeMatrix/Factorize
end

[mm,nn] = size(m_ts_jac);
if (length(inds) < nn & driver.oem.dofit)
  fprintf(1,'  length(inds) = %4i \n',length(inds));
  fprintf(1,'  size jac mm,nn = %4i %4i\n',mm,nn);
  disp('More jacobians than channels! will subset below')
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
wah = driver.rateset.unc_rates; [mgah,ngah] = size(wah);
if mgah == 1 | ngah == 1
  e0 = diag(driver.rateset.unc_rates(inds));
else
  e0 = driver.rateset.unc_rates(inds,inds);
end;  

% Error correlation matrix of observations (diagonal)
if mgah == 1 | ngah == 1
  i_e0_MatrOrArray = -1;     %% e0 = obs spectral uncertainty, is array
else
  i_e0_MatrOrArray = -1;     %% e0 = obs spectral uncertainty, is matrix
end

if i_e0_MatrOrArray < 0
  %% orig code
  se = e0 + fme;  
  se = se.*se;
  if isfield(aux,'all_obscov')
    %% this is in ../AIRS_new_random_scan_Aug2018/strow_override_defaults_latbins_AIRS.m
    se = aux.all_obscov;
  end
elseif i_e0_MatrOrArray > 0
  %% new code
  fme = diag(fme);
  fme = fme.*fme;
  if mgah == 1 | ngah == 1
    disp('  e0 = array ==> sent in an array of observational uncertainties')
    e0 = diag(e0);   %% sent in an array of uncertainties
    e0 = e0.*e0;
  else 
    disp('  e0 = matrix ==> sent in a matrix of observational uncertainties')
    e0 = e0;  %% sent cov matrix of obs uncertainties
  end
    se = e0 + fme;
end

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

%{
 inv_se = inv(se);      x0 = norm(eye(size(se)) - inv_se * se,'fro');
 inv_se = pinv(se);     x1 = norm(eye(size(se)) - inv_se * se,'fro');
 %inv_se = invillco(se); x2 = norm(eye(size(se)) - inv_se * se,'fro');
 inv_se = inverse(se);  x3 = norm(eye(size(se)) - inv_se * se,'fro');
 %fprintf(1,'norms = %8.6f %8.6f %8.6f %8.6f \n',[x0 x1 x2 x3]);
 fprintf(1,'norms = %8.6f %8.6f %8.6f \n',[x0 x1    x3]);
 error('kjsf')
%}

% Do this once to save time, assume diagonal, no need for pinv   ORIG 201
if invtype == 0
  inv_se = inv(se);         disp(' <<<<<<<< inv_se = inv(se)');            %%% NEW  pre Dec 2012
elseif invtype == 1
  inv_se = pinv(se);        disp(' <<<<<<<< inv_se = pinv(se)');           %%% NEW  post Dec 2012
elseif invtype == 2
  inv_se = invillco(se);    disp(' <<<<<<<< inv_se = invillco(se)');       %%% NEW  post Apr 2019
elseif invtype == 3
  inv_seF = factorize(se);  disp(' <<<<<<<< inv_se = factorize(se)');     %%% NEW  post Apr 2019
  inv_se = inverse(se);     disp(' <<<<<<<< inv_se = inverse(se)');       %%% NEW  post Apr 2019
  inv_se = inv_se * eye(size(inv_se));
end
%if invtype <= 2
  oo = find(isinf(inv_se) | isnan(inv_se)); inv_se(oo) = 0;
%end

% driver.oem.cov needs to be inverted, since it is literally covariance; 
% the tikonov matrices below ("s") need not be inverted
if invtype == 0
  rcov = inv(driver.oem.cov);
elseif invtype == 1
  rcov = pinv(driver.oem.cov);
elseif invtype == 2
  rcov = invillco(driver.oem.cov);
elseif invtype == 3
  rcovF = factorize(driver.oem.cov);
  rcov  = inverse(driver.oem.cov);
  rcov  = rcov * eye(size(rcov));
end

% Use following line for only Tikhonov reg THIS SHOULD NOT BE INVERTED
l = get_l(driver.jacobian.numlays,1);    
s = transpose(l)*l;

%% now build the Tikhonov regularization block matrix, using "s"
lenS = length(driver.jacobian.scalar_i);
%% default : only column/stemp jacs, layer T, layer WV
rc = blkdiag(zeros(lenS,lenS),driver.oem.alpha_water*s,driver.oem.alpha_temp*s);
if isfield(driver.oem,'alpha_ozone')
  rc = blkdiag(zeros(lenS,lenS),driver.oem.alpha_water*s,driver.oem.alpha_temp*s,driver.oem.alpha_ozone*s);
end

% if invtype == 3, rcov could be a class rather than double
switch driver.oem.reg_type
  case 'reg_and_cov'
    r = rcov + rc;
  case 'reg'
    r = rc;
  case 'cov'
    r = rcov;
  otherwise
    disp('Incorrect choice driver.oem.reg_type')
end

%whos r k inv_se
for ii = 1 : driver.oem.nloop
  % Do the retrieval inversion
  if invtype <= 2
    dx1 = r + k' * inv_se * k;
  elseif invtype == 3
    dx1 = inv_seF\k;
    dx1 = r + k'*dx1;
  end

  if invtype == 0
    dx1  = inv(dx1);
  elseif invtype == 1
    dx1  = pinv(dx1);
  elseif invtype == 2
    dx1  = invillco(dx1);
  elseif invtype == 3
    dx1 = double(dx1);
    dx1F = factorize(dx1);
    dx1  = inverse(dx1);
    dx1  = dx1 * eye(size(dx1));
  end
  if invtype <= 2
    dx2 = k' * inv_se * deltan - r*(xn-xb);
  elseif invtype == 3
    dx2 = k'/inv_seF * deltan - r*(xn-xb);
  end
  deltax = dx1*dx2; 

%{
  figure(1); plot(1:length(deltan),k);
  figure(2); pcolor(inv_se); shading flat; colorbar
  figure(2); plot(1:length(deltan),1./sqrt(diag(inv_se)),'b',1:length(deltan),-1./sqrt(diag(inv_se)),'b',1:length(deltan),deltan,'r'); grid
  figure(3); plot(dx2); title('dx2');
  figure(4); plot(dx1*dx2); title('DX1 * dx2');
  %error('lks')
%}

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
    chisqr(ii) = nansum(deltan'.*deltan');
   
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

%lala1 = diag(inv(AKstuff)*AKstuff);
%lala2 = diag(pinv(AKstuff)*AKstuff);
%lala3 = diag(invillco(AKstuff)*AKstuff);
%lala4 = diag(inv(U)*inv(L) * AKstuff);
%figure(5); plot(1:length(lala1),lala1,1:length(lala1),lala2,1:length(lala1),lala3); error('klfkf')

if invtype == 0
  errorx = inv(k' * inv_se * k + r);     %% decided pinv is too unstable, Aug 2018
elseif invtype == 1
  errorx = pinv(k' * inv_se * k + r);    %% decided pinv is too unstable, Nov 2013, but not much difference
elseif invtype == 2
  errorx = invillco(k' * inv_se * k + r);%% decided pinv is too unstable, Nov 2013, but not much difference
elseif invtype == 3
  errorx = double(k' * inv_se * k + r);  %% decided pinv is too unstable, May 2019, but not much difference
  errorx = inverse(errorx);
  errorx = errorx * eye(size(errorx));
elseif invtype == 9999
  AKstuff = (k' * inv_se * k + r);
  [L,U] = lu(AKstuff);
  errorx = inv(U)*inv(L);                %% trying LU
end

dofsx  = errorx * r; 
dofsx  = eye(size(dofsx)) - dofsx; 
dofs   = trace(dofsx);
cdofs  = diag(dofsx);                 %% so we can do cumulative d.of.f

% Gain is relative weight of first guess and observations
r_water  = r(driver.jacobian.water_i,driver.jacobian.water_i); 
r_temp   = r(driver.jacobian.temp_i,driver.jacobian.temp_i); 
if invtype == 0
  inv_r       = inv(r);
  inv_r_water = inv(r_water); 
  inv_r_temp  = inv(r_temp); 
elseif invtype == 1
  inv_r       = pinv(r);
  inv_r_water = pinv(r_water); 
  inv_r_temp  = pinv(r_temp); 
elseif invtype == 2
  inv_r       = invillco(r);
  inv_r_water = invillco(r_water); 
  inv_r_temp  = invillco(r_temp); 
elseif invtype == 3
  inv_r       = inverse(r);
  inv_r_water = inverse(r_water); 
  inv_r_temp  = inverse(r_temp); 
end

% inv operator seems OK for this matrix; if problems go back to pinv
k_water    = k(:,driver.jacobian.water_i); 
k_temp     = k(:,driver.jacobian.temp_i); 
if invtype == 0
  gain       = inv_r *k' * inv(k * inv_r * k' + se);
  gain_water = inv_r_water*k_water'*inv(k_water*inv_r_water*k_water'+se); 
  gain_temp  = inv_r_temp*k_temp'*inv(k_temp*inv_r_temp*k_temp'+se); 
elseif invtype == 1
  gain       = inv_r *k' * pinv(k * inv_r * k' + se);
  gain_water = inv_r_water*k_water'*pinv(k_water*inv_r_water*k_water'+se); 
  gain_temp  = inv_r_temp*k_temp'*pinv(k_temp*inv_r_temp*k_temp'+se); 
elseif invtype == 2
  gain       = inv_r *k' * invillco(k * inv_r * k' + se);
  gain_water = inv_r_water*k_water'*invillco(k_water*inv_r_water*k_water'+se); 
  gain_temp  = inv_r_temp*k_temp'*invillco(k_temp*inv_r_temp*k_temp'+se); 
elseif invtype == 3
  gain       = inv_r *k' * inverse(k * inv_r * k' + se);
  gain_water = inv_r_water*k_water'*inverse(k_water*inv_r_water*k_water'+se); 
  gain_temp  = inv_r_temp*k_temp'*inverse(k_temp*inv_r_temp*k_temp'+se); 
end

% Compute averaging kernel
ak = gain * k;   
ak_water = gain_water*k_water; 
ak_temp  = gain_temp*k_temp; 

if isfield(driver.oem,'alpha_ozone')
  r_ozone  = r(driver.jacobian.ozone_i,driver.jacobian.ozone_i); 
  k_ozone    = k(:,driver.jacobian.ozone_i);
  if invtype == 0 
    inv_r_ozone = inv(r_ozone);
  elseif invtype == 1 
    inv_r_ozone = pinv(r_ozone);
  elseif invtype == 2
    inv_r_ozone = invillco(r_ozone);
  elseif invtype == 3
    inv_r_ozone = inverse(r_ozone);
  end
  gain_ozone = inv_r_ozone*k_ozone'*inv(k_ozone*inv_r_ozone*k_ozone'+se);   
  ak_ozone = gain_ozone*k_ozone;
else
  ak_ozone = zeros(size(ak_water));
end
