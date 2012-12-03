function [rodgers_rate,errorx,dofs,gain,ak,inds,Se_errors] = rodgers(driver,aux_stuff)

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

addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/CLOUD
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
% Observation errors
ncerrors = aux_stuff.ncerrors * driver.oem.adjust_spectral_errorbars;

% Covariance (uncertainties/correlations) of measurements
lenr = length(inds);
fme  = ones(1,lenr)*driver.sarta_error;
fme  = diag(fme);          

Se_errors.ncerrors = ncerrors;
Se_errors.fmerrors = ones(size(ncerrors)) * driver.sarta_error;

% Obs errors with serial correlation correction
xe0 = diag(ncerrors(inds));        
if ~strcmp('junk',driver.oem.spectralcov_filename)
  spectral_cov = load(driver.oem.spectralcov_filename);
  dacov = spectral_cov.dacov;
  %% make sure there are some entries here!!!! ie not all zeros
  if nansum(diag(abs(dacov))) > 0
    e0 = dacov(inds,inds);
    %% now need to normalize
    xdiag = diag(xe0);
    xij = xdiag * xdiag';
    thediag = diag(e0);
    oo = find(isnan(e0(:))); e0(oo) = 0;
    oo = find(isinf(e0(:))); e0(oo) = 0;
    whos e0 xe0 xdiag
    e0 = e0 * nanmean(xdiag)/nanmean(diag(e0));
  else
    e0 = xe0;
    fprintf(1,'Warning : not using spectral covs from %s as they are NaNs \n',driver.oem.spectralcov_filename)
  end
else  %% just use ncerrors
  e0 = xe0;
end

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
oo = find(isinf(inv_se) | isnan(inv_se)); inv_se(oo) = 0;

% Apply regularization multiplier
if driver.oem.diag_only
  % Paul's suggestion to access the diagonal
  r(1:size(r,1)+1:end) = r(1:size(r,1)+1:end) * driver.oem.lambda;

else
  %%%%%%%%%%%%%%%%%%%% this is for QST ie the SINGLE COLUMN gases and SST
  % Multiply entire matrix, but in blocks of column(QST), qlays(Q1 .. QN), templays
  if length(driver.oem.lambda_qst) == 1
    %% multiply all R(i,j) of the qst part of the matrix, by one number
    r(driver.jacobian.iqst,driver.jacobian.iqst)     = r(driver.jacobian.iqst,driver.jacobian.iqst) * driver.oem.lambda_qst;
  elseif length(driver.oem.lambda_qst) == length(driver.jacobian.iqst)
    %% multiply all R(i,i) of the qst part of the matrix, ,* specified diagnol, making this diagnol
    r(driver.jacobian.iqst,driver.jacobian.iqst)     = r(driver.jacobian.iqst,driver.jacobian.iqst) .* diag(driver.oem.lambda_qst);
  elseif length(driver.oem.lambda_qst) > length(driver.jacobian.iqst)
    %% matrix multiply all R(i,j) of the qst part of the matrix, by specified diagnol
    junk = driver.oem.lambda_qst(1:length(driver.jacobian.iqst));
    junk = diag(junk);
    r(driver.jacobian.iqst,driver.jacobian.iqst)     = r(driver.jacobian.iqst,driver.jacobian.iqst) * junk;
  end

  %%%%%%%%%%%%%%%%%%%% this is for Q(z) ie the MULTICOLUMN gases
  for ii = 1 : driver.jacobian.numQlays
    junk = ['lala = driver.oem.lambda_Q' num2str(ii) ';']; eval(junk);
    junk = ['mama = driver.jacobian.iQ' num2str(ii) ';'];  eval(junk);
    if length(lala) == 1
      %% multiply all R(i,j) of the iQ part of the matrix, by one number
      r(mama,mama) = r(mama,mama) * lala;
    elseif length(lala) == length(mama)
      %% multiply all R(i,i) of the iQ part of the matrix, .* specified diagnol, making this diagnol
      r(mama,mama) = r(mama,mama)  .* diag(lala);
    elseif length(lala) > length(mama)
      %% multiply all R(i,j) of the iQ part of the matrix, by specified diagnol
      junk = lala(1:length(mama));
      junk = diag(junk);
      r(mama,mama) = r(mama,mama)  * (junk);
    end
  end

  %%%%%%%%%%%%%%%%%%%% this is for T(z) ie the MULTICOLUMN temperature
  lala = driver.oem.lambda_temp;
  mama = driver.jacobian.itemp;
  if length(driver.oem.lambda_temp) == 1
    %% multiply all R(i,j) of the T part of the matrix, by one number
    r(driver.jacobian.itemp,driver.jacobian.itemp)   = r(driver.jacobian.itemp,driver.jacobian.itemp) * driver.oem.lambda_temp;
  elseif length(lala) == length(mama)
    %% multiply all R(i,i) of the T part of the matrix, .* specified diagnol, making this diagnol
    r(driver.jacobian.itemp,driver.jacobian.itemp)   = r(driver.jacobian.itemp,driver.jacobian.itemp) .* diag(driver.oem.lambda_temp);
  elseif length(lala) > length(mama)
    %% multiply all R(i,j) of the iQ part of the matrix, by specified diagnol
    junk = lala(1:length(mama));
    junk = diag(junk);
    r(mama,mama) = r(mama,mama)  * (junk);
  end

  % Now add diagonal, but only for T and water
  n = length(driver.jacobian.iQ1);
  m = length(driver.jacobian.itemp);
  ti = driver.jacobian.itemp;
  for i=1:m
    r(ti(i),ti(i)) = r(ti(i),ti(i)) + driver.oem.lambda;
  end

  for ii = 1 : driver.jacobian.numQlays
    junk = ['qi = driver.jacobian.iQ' num2str(ii) ';']; eval(junk);
    for i=1:n
      r(qi(i),qi(i)) = r(qi(i),qi(i)) + driver.oem.lambda;
    end
  end

  %% fine add it on for QST as well (Aug 2012)
  n = length(driver.jacobian.qstYesOrNo);
  for i=1:n
    r(i,i) = r(i,i) + driver.oem.lambda;
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

%keyboard
