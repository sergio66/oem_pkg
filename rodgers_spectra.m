function [rodgers_rate,errorx,dofs,cdofs,gain,ak,inds,r,se,inv_se,se_errors] = rodgers(driver,aux_stuff)

%---------------------------------------------------------------------------
% OEM retrieval for REGULAR SPECTRA so need klayers and sarta exec
%   so y(x) = klayers/sarta, to compare to yIN
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
%   se                 = actual channel covariance matrix used for spectral obs (1x2378 or 1x8461)
%   inv_se             = actual channel inverse covariance matrix used
%   se_errors          = actual channel uncertainties used
%
%---------------------------------------------------------------------------
addpath /asl/matlib/aslutil/
rtpop = mktemp('junk.op.rtp');
rtprp = mktemp('junk.rp.rtp');

% Load and subset relaxation matrix fpr PARAMETERS
r = load(driver.oem.cov_filename,'cov');
r = r.cov(driver.jacindex,driver.jacindex);

% Jacobians
m_ts_jac = aux_stuff.m_ts_jac; 
[junkMM,junkNN] = size(m_ts_jac);

% this is only for checking sizes of the jacs, to see what "instrument" they correspond to
if junkMM ~= 2378 & junkMM ~= 8461
  error('can only handle AIRS or IASI')
end
clear junkMM junkNN

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

se_errors.ncerrors = ncerrors;
se_errors.fmerrors = ones(size(ncerrors)) * driver.oem.sarta_error;

%% get 2378x2378 spectral cov matrix
e0 = get_spectral_covariances(driver,ncerrors,inds);

% Error correlation matrix of observations (diagonal)
se = e0 + fme;  
se = se.*se;

% xb is the apriori
[zz1,zz2] = size(xb);
% Linearization point = zero, assuming fits are linear
xn = xb;              %% after July 2013

% Form k matrix (Jacobians)
k = m_ts_jac(inds,:);
[mm,nn] = size(k);

% Form y - F(xa)
tempx = run_klayers_sarta(driver,xn,rtpop,rtprp,0,-1);
deltaY = driver.rateset.rates(inds) - tempx(inds);
deltaY0 = deltaY;
plot(1:length(deltaY),deltaY,1:length(deltaY),deltaY,'r');
plot(driver.f(inds),deltaY,'r'); pause(0.1)

% Do this once to save time, assume diagonal, no need for pinv   ORIG 201
% ???????????????
% inv_se = diag(1./diag(se)); disp(' >>>>>>>> inv_se = diag(1./diag(se))');  %%% ORIG pre Dec 012
inv_se = pinv(se);          disp(' <<<<<<<< inv_se = pinv(se)');           %%% NEW  post Dec 2012
oo = find(isinf(inv_se) | isnan(inv_se)); inv_se(oo) = 0;
% ???????????????

if driver.oem.regularizationVScovariances == 'R' | driver.oem.regularizationVScovariances == 'r'
  % Apply (smoothing) regularization multiplier 
  pm1 = r;
  pm1 = r*r';  % AVT do we need this? This makes pm1 the square of r.
  r = regularization_multiplier(pm1,driver);
  % override_cov_r       % do we want to couple SST and TS,TT or CO2 and TS,TT???
  % override_cov_rVERS2  % do we want to use COV from ERA?
elseif driver.oem.regularizationVScovariances == 'C' | driver.oem.regularizationVScovariances == 'c'
  % Apply covariance matrix, which has correlations
  r = geophysical_covariance(driver);
else
  error('need driver.oem.regularizationVScovariances == r or R or c or C')
end

for ii = 1 : driver.oem.nloop
  % Do the retrieval inversion
  dx1    = r + k' * inv_se * k; 
  dx1    = pinv(dx1);    

  dx2    = k' * inv_se * deltaY - r*(xn-xb);
  deltax = dx1*dx2; 

  % Update first guess with deltax changes
  rodgers_rate = real(xn + deltax);  
  xn = rodgers_rate;

  xsave(ii,:) = rodgers_rate;

  if ii < driver.oem.nloop
    channel_errors = ncerrors(inds);
    deltaYIN = deltaY;

    % Form the computed rates; also see lines 108-121 of oem_lls.m
    tempx = run_klayers_sarta(driver,xn,rtpop,rtprp,ii,-1);
    %[tempx,jacnew] = run_klayers_sarta(driver,xn,rtpop,rtprp,ii,+1);
    %jacnew = jacnew(inds,:);
    %k = jacnew;
 
    % Compute chisqr, and new deltaY
    deltaY = driver.rateset.rates - tempx;
    deltaY = deltaY(inds);
    chisqr(ii) = sum(deltaY.*deltaY);

    figure(2); clf
    plot(1:length(deltaY),deltaYIN,1:length(deltaY),deltaY,'r'); 
    plot(driver.f(inds),deltaYIN,driver.f(inds),deltaY,'r',driver.f(inds),deltaY0,'g--',...
         driver.f(inds),channel_errors,'k',driver.f(inds),-channel_errors,'k'); 
    title(['obs - fit at iteration ' num2str(ii)]); pause(0.1)
  end

end

best = find(chisqr ==  min(chisqr),1);
rodgers_rate = xsave(best,:);

if driver.oem.nloop > 1
  disp('printing out successive chisqr values (upto N-1 th iterate) ...')
  chisqr
end

% Error analysis and diagnostics
errorx = pinv(k' * inv_se * k + r);
dofsx  = errorx * r; 
dofsx  = eye(size(dofsx)) - dofsx; 
dofs   = trace(dofsx);
cdofs  = diag(dofsx);

% Gain is relative weight of first guess and observations
inv_r = pinv(r);

% inv operator seems OK for this matrix; if problems go back to pinv
gain  = inv_r *k' * inv(k * inv_r * k' + se);

% Compute averaging kernel
ak = gain * k;   

rmer = ['!/bin/rm ' rtpop ' ' rtprp];
eval(rmer);