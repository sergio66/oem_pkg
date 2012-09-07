function [rodgers_rate,errorx,dofs,gain,ak,inds] = rodgers(driver,aux_stuff)

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
%---------------------------------------------------------------------------

% Load and subset relaxation matrix
r = load(driver.oem.cov_filename,'cov');
r = r.cov(driver.jacindex,driver.jacindex);

% Jacobians
m_ts_jac = aux_stuff.m_ts_jac; [junkMM,junkNN] = size(m_ts_jac);

addpath /home/sergio/MATLABCODE
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

%o3chans = find(fairs > 1000 & fairs < 1080); o3chans  = intersect(o3chans,gg);
%cfcchans = find(fairs > 843 & fairs < 848);  cfcchans = intersect(cfcchans,gg);
%inds = [inds o3chans cfcchans];

fprintf(1,'orig chan length = %4i \n',length(inds));
cfcchans = find(fairs >= 842.5 & fairs < 849.0);  cfcchans = intersect(cfcchans,gg); 
fprintf(1,'added on %4i CFC chans \n',length(cfcchans));
inds = [inds cfcchans];

% Apriori state
xb       = aux_stuff.xb;
% Observation errors
ncerrors = aux_stuff.ncerrors;

% Covariance (uncertainties/correlations) of measurements
lenr = length(inds);
fme  = ones(1,lenr)*driver.sarta_error;
% ssm xjunk = ncerrors(inds);        
% ssm fme(o3chans) = 1e10*xjunk(o3chans);

fme  = diag(fme);          

%o3chans = find(fairs(inds) > 1000 & fairs(inds) < 1080);
%cfcchans = find(fairs(inds) > 843 & fairs(inds) < 848);

% Obs errors with serial correlation correction
e0 = diag(ncerrors(inds));        

% Error correlation matrix of observations (diagonal)
se = e0 + fme;  
se = se.*se;

% xb is the apriori
[l1,l2] = size(xb);
% Linearization point = zero, assuming fits are linear
xn = zeros(l1,l2);

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

% Apply regularization multiplier
if driver.oem.diag_only
  % Paul's suggestion to access the diagonal
  r(1:size(r,1)+1:end) = r(1:size(r,1)+1:end) * driver.oem.lambda;
else
  % Multiply entire matrix, but in blocks of gas, water, temp
  r(driver.jacobian.igas,driver.jacobian.igas)     = r(driver.jacobian.igas,driver.jacobian.igas) * driver.oem.lambda_gas;
  % ssm r(2,2) = r(2,2) + 100/driver.oem.lambda;  %% <<<<<<<<<<<------ new
  % ssm r(1:5,1:5)
  r(driver.jacobian.iwater,driver.jacobian.iwater) = r(driver.jacobian.iwater,driver.jacobian.iwater) * driver.oem.lambda_water;
  r(driver.jacobian.itemp,driver.jacobian.itemp)   = r(driver.jacobian.itemp,driver.jacobian.itemp) * driver.oem.lambda_temp;

  % When only a single lambda
  % r = r * driver.oem.lambda;

  % Now add diagonal, but only for T and water
  n = length(driver.jacobian.iwater);
  m = length(driver.jacobian.itemp);
  wi = driver.jacobian.iwater;
  ti = driver.jacobian.itemp;
  for i=1:n
    r(wi(i),wi(i)) = r(wi(i),wi(i)) + driver.oem.lambda;
  end
  for i=1:m
    r(ti(i),ti(i)) = r(ti(i),ti(i)) + driver.oem.lambda;
  end
end

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
