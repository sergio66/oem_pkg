function [rodgers_rate,errorx,dofs,gain,ak,inds,Se_errors,r,inv_se] = rodgers(driver,aux_stuff)

%---------------------------------------------------------------------------
% OEM retrieval
%---------------------------------------------------------------------------
% Notation consistent with 
% Tilman Steck, Methods for determining regularization for atmospheric 
%               retieval problems, Appl Optics, v41, pg 1788 (2002)
%
% jacobian           = k
% obs error matrix   = se
% param error matrix = r0,r
% apriori            = xset
% inital value       = xn
% final value        = xnp1 = xn + deltax
% deg of freedom     = dofs
% gain matrix        = gain
% averaging kernel   = ak
% inds               = what channels used
% Se_errors          = uncertainties used in building up Se matrix (from obs and forward model)
%---------------------------------------------------------------------------

% Load and subset relaxation matrix
r = load(driver.oem.cov_filename,'cov');
r = r.cov(driver.jacindex,driver.jacindex);
rInit = r;

% Jacobians
m_ts_jac = aux_stuff.m_ts_jac; [junkMM,junkNN] = size(m_ts_jac);

%addpath /home/sergio/MATLABCODE
%addpath /home/sergio/MATLABCODE/CLOUD
if junkMM == 2378
  fairs = instr_chans('airs');
  gg    = dogoodchan;
elseif junkMM == 8461
  fairs = instr_chans('iasi');
  gg    = 1:8461;
else
  error('can only handle AIRS or IASI')
end

% Index of frequencies used
inds     = driver.jacobian.chanset;

% Apriori state
xb       = aux_stuff.xb;
% Observation errors : can adjust them with a scalar 
ncerrors = aux_stuff.ncerrors * driver.oem.adjust_spectral_errorbars;

% Covariance (uncertainties/correlations) of measurements
lenr = length(inds);
fme  = ones(1,lenr)*driver.oem.sarta_error;
fme  = diag(fme);          

Se_errors.ncerrors = ncerrors;
Se_errors.fmerrors = ones(size(ncerrors)) * driver.oem.sarta_error;

%% get 2378x2378 spectral cov matrix
e0 = get_spectral_covariances(driver,ncerrors,inds);

% Error correlation matrix of observations (diagonal)
se = e0 + fme;  
se = se.*se;

% xb is the apriori
[zz1,zz2] = size(xb);
% Linearization point = zero, assuming fits are linear
xn = zeros(zz1,zz2);

% Form k matrix (Jacobians)
k = m_ts_jac(inds,:);
[mm,nn] = size(k);

% Form y - F(xa)
fx = zeros(size(driver.rateset.rates));
for iy = 1 : length(xn)
   fx = fx + (xn(iy)*m_ts_jac(:,iy));
end
deltan = driver.rateset.rates(inds) - fx(inds);

% Do this once to save time, assume diagonal, no need for pinv
inv_se = diag(1./diag(se));
oo = find(isinf(inv_se) | isnan(inv_se)); inv_se(oo) = 0;

% Apply regularization multiplier 
r = regularization_multiplier(r,driver);
% override_cov_r       % do we want to couple SST and TS,TT or CO2 and TS,TT???
% override_cov_rVERS2  % do we want to use COV from ERA?

% Do the retrieval inversion
dx1    = r + k' * inv_se * k; 
dx1    = pinv(dx1);    
dx2    = k' * inv_se * deltan - r*(xn-xb);
deltax = dx1*dx2; 

% Update first guess with deltax changes
rodgers_rate = real(xn + deltax)';  

% Error analysis and diagnostics
errorx = pinv(k' * inv_se * k + r);
dofsx  = errorx * r; 
dofsx  = eye(size(dofsx)) - dofsx; 
dofs   = trace(dofsx);

% Gain is relative weight of first guess and observations
inv_r = pinv(r);

% inv operator seems OK for this matrix; if problems go back to pinv
gain  = inv_r *k' * inv(k * inv_r * k' + se);

% Compute averaging kernel
ak = gain * k;   

