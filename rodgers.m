function [rodgers_rate,errorx,dofs,cdofs,gain,ak,r,se,inv_se,se_errors,ak_water,ak_temp,ak_ozone,bestloop,deltan00] = rodgers(driver,aux)

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
%   bestloop           = iteration number where minimum chisqr was found; the oem params saved at this point
%
%---------------------------------------------------------------------------
% Jacobians
m_ts_jac = aux.m_ts_jac;

% Index of frequencies used
inds     = driver.jacobian.chanset;
%inds     = inds(inds <= 600);  'to only do 15 um chans in rodgers.m'

invtype = 0;  %% inv
invtype = 1;  %% pinv     BEST *****
invtype = 2;  %% S. Rump      invillco   addpath /home/sergio/MATLABCODE/IntLab
invtype = 3;  %% T. A. Davis  factorize  addpath /home/sergio/MATLABCODE/FactorizeMatrix/Factorize
invtype = 4;  %% for Se do a ridge regression    Se --> Senew = Se + delta I
invtype = 5;  %% for Se do a minimum eigenvalue  Se --> Senew = Se + blah (eig > minimum)
if ~isfield(aux,'invtype')
  aux.invtype = 1;   %% default is to use pinv
end
invtype = aux.invtype;
if invtype < 0 | invtype > 5
  error('need invtype between 0 and 5')
end
fprintf(1,'inverse of matrices using method (0) inv (1) PINV (default) (2) invillco (3) factorize (4) Se RR (5) Se ME : %2i \n',invtype);

% max condition number for invtype == 4
kmax = 1000;  
kmax = 10000;  
kmax = 100000;  
kmax = 1e4;  %% works pretty good  
kmax = 1e1; 

% min eigenvalue for invtype == 5
sigmin = 1.0e-12;
sigmin = 1.0e-10; %% works pretty good
sigmin = 1.0e-16;

addpath /home/sergio/MATLABCODE
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

% Apriori state; make sure it has been correctly normalized before being used here!
xb       = aux.xb;

% Covariance (uncertainties/correlations) of measurements
lenr = length(inds);
fme  = ones(1,lenr)*driver.oem.sarta_error;
fme  = diag(fme);          

sizer = size(driver.rateset.rates);
se_errors.fmerrors = ones(sizer) * driver.oem.sarta_error;

%% get 2378x2378 spectral cov matrix
wah = driver.rateset.unc_rates; [mgah,ngah] = size(wah);
%size(driver.rateset.unc_rates)
%size(inds)
if mgah == 1 | ngah == 1
  e0 = diag(driver.rateset.unc_rates(inds));
else
  e0 = driver.rateset.unc_rates(inds,inds);
end;  

% Error correlation matrix of observations (diagonal)
if mgah == 1 | ngah == 1
  i_e0_MatrOrArray = -1;     %% e0 = obs spectral uncertainty, is array
else
  i_e0_MatrOrArray = +1;     %% e0 = obs spectral uncertainty, is matrix
end

if i_e0_MatrOrArray < 0
  %% orig code, send in vector of spectral uncertainty ... so turn it into matrix
  se = e0 + fme;  
  se = se.*se;
  if isfield(aux,'all_obscov')
    %% this is in ../AIRS_new_random_scan_Aug2018/strow_override_defaults_latbins_AIRS.m
    disp('using aux.all_obscov')
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
xnIN = xn;

% Form k matrix (Jacobians)
k = m_ts_jac(inds,:);
[mm,nn] = size(k);

iAddXB = -1; %% new, does this really makes more sense see eg anomaly_0dayavg_resultsXloop3try2?????
iAddXB = +1; %% orig, wierd but I think it is ok as you need delta0 = obs - tracegas_offset = obs-f(x0) = obs - f(xb)
if iAddXB > 0
  nyuk = find(abs(xn) > eps);
  %[nyuk xn(nyuk)]
  % Form y - F(xa), this is orig code but a little wierd!!!!!!
  tracegas_offset = zeros(size(driver.rateset.rates));
  for iy = 1 : length(xn)
     tracegas_offset = tracegas_offset + (xn(iy)*m_ts_jac(:,iy));
  end
  iJUNK = [driver.jacobian.scalar_i  driver.jacobian.water_i(1) driver.jacobian.temp_i(1) driver.jacobian.ozone_i(1)];
  disp('   scalar/WV/T/O3 xb');
  disp('   xb    |qrenorm   |  xb.*qrenorm')
  disp('------------------------------------')
  fprintf(1,'%8.4f | %8.4f | %8.4f \n', [xn(iJUNK)   driver.qrenorm(iJUNK)'  xn(iJUNK).*driver.qrenorm(iJUNK)']')
  disp('------------------------------------')

  tracegas_offset00 = tracegas_offset;
  deltan00 = driver.rateset.rates - tracegas_offset00;    %%% << this is what we are fitting, all chans >>
  deltan   = deltan00(inds);                 %%% << this is what we are fitting, selected chans >>
  deltan0  = deltan;
else
  tracegas_offset = zeros(size(driver.rateset.rates));
  tracegas_offset00 = tracegas_offset;
  deltan00 = driver.rateset.rates - tracegas_offset00;    %%% << this is what we are fitting, all chans >>
  deltan   = deltan00(inds);                 %%% << this is what we are fitting, selected chans >>
  deltan0  = deltan;
end

% hdffile = '/home/sergio/MATLABCODE/airs_l1c_srf_tables_lls_20181205.hdf';   % what he gave in Dec 2018
% vchan2834 = hdfread(hdffile,'freq');
% f = vchan2834;
% load sarta_chans_for_l1c.mat
% f = f(ichan);
% f = f(inds);
if length(driver.rateset.rates) == 2645
  f = instr_chans2645;
elseif length(driver.rateset.rates) == 2378
  f = instr_chans;
elseif length(driver.rateset.rates) == 1305
  f = instr_chans('cris1305');
else
  error('oooorrr is this AIRS 2378 or 2465 or Cris 1305?')
end

figure(12); plot(1:length(tracegas_offset),tracegas_offset,1:length(tracegas_offset),m_ts_jac(:,1:3)); hl = legend('tracegas offset','CO2 jac','N2O jac','CH4 jac','location','best','fontsize',10);

figure(1); plot(f(inds),driver.rateset.rates(inds),'b.-',f(inds),tracegas_offset(inds),'g.-',f(inds),driver.rateset.rates(inds) - tracegas_offset00(inds),'c.-','linewidth',2); 
  plotaxis2; title('in oem\_pkg/rodgers.m : nyuk'); 
  hl = legend('input rates','trace gas jacs offset','signal''= to fit b-g','location','best');

[mmm,nnn] = size(m_ts_jac);
if nnn == 66
  wahCO2_ST = m_ts_jac(inds,[1 6]);
  wahWV = m_ts_jac(inds,(1:20)+6+00);
  wahT  = m_ts_jac(inds,(1:20)+6+20);
  wahO3 = m_ts_jac(inds,(1:20)+6+40);
  figure(1); plot(f(inds),driver.rateset.rates(inds) - tracegas_offset00(inds),'kx-',f(inds),sum(wahWV'),f(inds),sum(wahT'),f(inds),sum(wahO3'),f(inds),wahCO2_ST(:,1),f(inds),wahCO2_ST(:,2),'linewidth',2); 
    plotaxis2; title('in oem\_pkg/rodgers.m : nyuk'); 
    hl = legend('signal''= to fit b-g','WVjac','Tjac','O3jac','CO2jac','STjac','location','best','fontsize',10);
  figure(1); plot(f(inds),sum(wahWV'),f(inds),sum(wahT'),f(inds),sum(wahO3'),f(inds),wahCO2_ST(:,1),f(inds),wahCO2_ST(:,2),'linewidth',2); 
    plotaxis2; title('in oem\_pkg/rodgers.m : nyuk'); 
    hl = legend('WVjac','Tjac','O3jac','CO2jac','STjac','location','best','fontsize',10);

end
  
if iAddXB > 0
  %indsy791 = find(f >= 790,1); indsy791 = sort([inds; (indsy791-25:indsy791+25)']);

  %driver.oem.doplots = 1  
  if driver.oem.doplots > 0
    iTRPorSTD = +49;
    iTRPorSTD = +1;

    figure(2); plot(f(inds),m_ts_jac(inds,1)); title('Should be CO2 jac')
    figure(2); plot(f(inds),m_ts_jac(inds,1),'.-'); title('Should be CO2 jac'); xlim([640 840]); grid; grid minor
  
    if iTRPorSTD == 49
      miaow = load('/asl/s1/sergio/AIRSPRODUCTS_JACOBIANS/STD/g2_jac.mat');
    elseif iTRPorSTD == 1
      miaow = load('/asl/s1/sergio/AIRSPRODUCTS_JACOBIANS/TRP/g2_jac.mat');
    end
    figure(2); plot(f(inds),m_ts_jac(inds,1),'.-',miaow.fout,sum(miaow.jout')*2.2/370); title('Should be CO2 jac'); xlim([640 840]); grid; grid minor
      hl = legend('input jac','from STD/g2\_jac.mat','location','best','fontsize',10);
  
    if iTRPorSTD == 49
      miaow = load('/asl/s1/sergio/AIRSPRODUCTS_JACOBIANS/STD/g6_jac.mat');
    elseif iTRPorSTD == 1
      miaow = load('/asl/s1/sergio/AIRSPRODUCTS_JACOBIANS/TRP/g6_jac.mat');
    end
    iCH4 = 3; %% anomaly tile spectra
    figure(3); plot(f(inds),m_ts_jac(inds,iCH4),'.-',miaow.fout,sum(miaow.jout')*5/1860); title('Should be CH4 jac'); xlim([640 1340]); grid; grid minor
      hl = legend('input jac','from STD/g6\_jac.mat','location','best','fontsize',10);

    if iTRPorSTD == 49
      miaow = load('/asl/s1/sergio/AIRSPRODUCTS_JACOBIANS/STD/surface_jac.mat');
    elseif iTRPorSTD == 1
      miaow = load('/asl/s1/sergio/AIRSPRODUCTS_JACOBIANS/TRP/surface_jac_new.mat');
    end
    iST = 4; %% anomaly tile spectra
    figure(4); plot(f(inds),m_ts_jac(inds,iST),'.-',miaow.fout,miaow.jsurface(:,1)*0.1); title('Should be ST jac'); xlim([640 1340]); grid; grid minor
      hl = legend('input jac','from STD/surface\_jac.mat','location','best','fontsize',10);
  
    if iTRPorSTD == 49
      miaow1   = load('/asl/s1/sergio/AIRSPRODUCTS_JACOBIANS/STD/g1_jac.mat');
      miaow101 = load('/asl/s1/sergio/AIRSPRODUCTS_JACOBIANS/STD/g101_jac.mat');
      miaow102 = load('/asl/s1/sergio/AIRSPRODUCTS_JACOBIANS/STD/g102_jac.mat');
      miaow = miaow1;
      miaow.jout = miaow1.jout + miaow101.jout + miaow102.jout;
    elseif iTRPorSTD == 1
      miaow1   = load('/asl/s1/sergio/AIRSPRODUCTS_JACOBIANS/TRP/g1_jac_new.mat');
      miaow101 = load('/asl/s1/sergio/AIRSPRODUCTS_JACOBIANS/TRP/g101_jac_new.mat');
      miaow102 = load('/asl/s1/sergio/AIRSPRODUCTS_JACOBIANS/TRP/g102_jac_new.mat');
      miaow103 = load('/asl/s1/sergio/AIRSPRODUCTS_JACOBIANS/TRP/g103_jac_new.mat');
      miaow = miaow1;
      miaow.jout = miaow1.jout + miaow101.jout + miaow102.jout + miaow103.jout;
    end
    iWV = (01:20)+4; %% anomaly tile spectra
    figure(5); plot(f(inds),sum(m_ts_jac(inds,iWV),2),'.-',miaow.fout,sum(miaow.jout')*0.01); title('Should be WV jac'); xlim([640 1340]); grid; grid minor
      hl = legend('input jac','from STD/g1\_jac.mat','location','best','fontsize',10);

    if iTRPorSTD == 49  
      miaow   = load('/asl/s1/sergio/AIRSPRODUCTS_JACOBIANS/STD/temp_jac.mat');
    elseif iTRPorSTD == 1  
      miaow   = load('/asl/s1/sergio/AIRSPRODUCTS_JACOBIANS/TRP/temp_jac_new.mat');
    end
    iTz = 20+iWV; %% anomaly tile spectra
    figure(6); plot(f(inds),sum(m_ts_jac(inds,iTz),2),'.-',miaow.fout,sum(miaow.jtemp')*0.01); title('Should be T jac'); xlim([640 1340]); grid; grid minor
      hl = legend('input jac','from STD/temp\_jac.mat','location','best','fontsize',10);
  
    if iTRPorSTD == 49  
      miaow   = load('/asl/s1/sergio/AIRSPRODUCTS_JACOBIANS/STD/g3_jac.mat');
    elseif iTRPorSTD == 1
      miaow   = load('/asl/s1/sergio/AIRSPRODUCTS_JACOBIANS/TRP/g3_jac.mat');
    end
    iO3 = 40+iWV; %% anomaly tile spectra
    figure(7); plot(f(inds),sum(m_ts_jac(inds,iO3),2),'.-',miaow.fout,sum(miaow.jout')*0.01); title('Should be O3 jac'); xlim([640 1340]); grid; grid minor
      hl = legend('input jac','from STD/g3\_jac.mat','location','best','fontsize',10);
  pause(0.1)
  end
end
%disp('nyuk rodgers.m ret to continue'); pause

%pause(0.1);

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
  disp(' <<<<<<<< inv_se = inv(se)');            %%% NEW  pre Dec 2012
  inv_se = inv(se);         
elseif invtype == 1
  disp(' <<<<<<<< inv_se = pinv(se) DEFAULT');   %%% NEW  post Dec 2012
  inv_se = pinv(se);        
elseif invtype == 2
  disp(' <<<<<<<< inv_se = invillco(se)');       %%% NEW  post Apr 2019
  inv_se = invillco(se);    
elseif invtype == 3
  disp(' <<<<<<<< inv_se = inverse(se)');        %%% NEW  post Apr 2019
  inv_seF = factorize(se);  
  inv_se = inverse(se);     
  inv_se = inv_se * eye(size(inv_se));
elseif invtype == 4
  disp(' <<<<<<<< inv_se = inverse_ridge_regression_matrix(se)');            %%% NEW  pre Dec 2012
  inv_se = inverse_ridge_regression_matrix(se,kmax);   
elseif invtype == 5
  disp(' <<<<<<<< inv_se = inverse_minimum_eigenvalue_matrix(se)');            %%% NEW  pre Dec 2012
  kmaxrange   = 2 : 1 : 12;   kmaxrange = 10.^kmaxrange;
  sigminrange = -16 : 1 : -8; sigminrange = 10.^sigminrange;

  kmaxrange   = 2 : 0.25 : 12;   kmaxrange = 10.^kmaxrange;
  sigminrange = -16 : 0.25 : -8; sigminrange = 10.^sigminrange;

  kmaxrange   = 2 : 0.25 : 12;   kmaxrange = 10.^kmaxrange;
  sigminrange = -20 : 1 : -12; sigminrange = 10.^sigminrange;

  inv_se = inverse_minimum_eigenvalue_matrix_optim(se,kmaxrange,sigminrange,'Se');   
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
elseif invtype == 4
  rcov = inverse_ridge_regression_matrix(driver.oem.cov,kmax);
elseif invtype == 5
  rcov = inverse_minimum_eigenvalue_matrix_optim(driver.oem.cov,kmaxrange,sigminrange,'driver.oem.cov');
end

% Use following line for only Tikhonov reg THIS SHOULD NOT BE INVERTED
l = get_l(driver.jacobian.numlays,1);    
s = transpose(l)*l;

%% now build the Tikhonov regularization block matrix, using "s"
lenS = length(driver.jacobian.scalar_i);
%% default : only column/stemp jacs, layer WV
rc = blkdiag(zeros(lenS,lenS),driver.oem.alpha_water*s);

%% always assumes you want to fit for WV ... may want to keep T fixed 9so no fit) and also may not want to fit O3
if isfield(driver.oem,'alpha_temp') & isfield(driver.oem,'alpha_ozone')
  rc = blkdiag(zeros(lenS,lenS),driver.oem.alpha_water*s,driver.oem.alpha_temp*s,driver.oem.alpha_ozone*s);
elseif isfield(driver.oem,'alpha_temp') & ~isfield(driver.oem,'alpha_ozone')
  rc = blkdiag(zeros(lenS,lenS),driver.oem.alpha_water*s,driver.oem.alpha_temp*s);
elseif ~isfield(driver.oem,'alpha_temp') & isfield(driver.oem,'alpha_ozone')
  rc = blkdiag(zeros(lenS,lenS),driver.oem.alpha_water*s,driver.oem.alpha_ozone*s);
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
chisqr0 = nansum(deltan'.*deltan');

%whos rcov rc r k inv_se
%disp('rodgers.m 1'); keyboard_nowindow

for ii = 1 : driver.oem.nloop
  % Do the retrieval inversion
  if invtype ~= 3
    dx1 = r + k' * inv_se * k;
    figure(7); imagesc(log10(abs(r))); title('r'); colorbar
    figure(8); imagesc(log10(abs(k))); title('k'); colorbar
    figure(9); imagesc(log10(abs(inv_se))); title('inv_se'); colorbar
    pause(0.1)

  elseif invtype == 3
    dx1 = inv_seF\k;
    dx1 = r + k'*dx1;
  end

  dx1_0 = dx1;
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
  elseif invtype == 4
    dx1  = inverse_ridge_regression_matrix(dx1,kmax);
  elseif invtype == 5
    dx1  = inverse_minimum_eigenvalue_matrix_optim(dx1,kmaxrange,sigminrange,'dx1');
  end
  if invtype ~= 3
    dx2 = k' * inv_se * deltan - r*(xn-xb);
  elseif invtype == 3
    dx2 = k'/inv_seF * deltan - r*(xn-xb);
  end
  deltax = dx1*dx2;
  figure(4); plot(diag(dx1)); colorbar;                              title('log10(dx1)');
  figure(5); imagesc(log10(abs(dx1))); colorbar;                     title('dx1');
  figure(6); plot(dx2);                                              title('dx2');
  figure(7); plot(f(inds),deltan); plotaxis2;                        title('deltaBT to fit')
  figure(8); plot(deltax.*driver.qrenorm'); plotaxis2; grid minor;   title('deltax.*qrenorm')

  iDebug = +1;
  iDebug = -1;
  if iDebug > 0

    %addpath /home/sergio/MATLABCODE; keyboard_nowindow
    figure(1); plot(f,k); grid
    figure(1); plot(f,k(:,1:5)); grid
%  figure(1); plot(f,k(:,1)); grid
    figure(2); pcolor(inv_se); shading flat; colorbar
    figure(2); plot(1:length(deltan),1./sqrt(diag(inv_se)),'b',1:length(deltan),-1./sqrt(diag(inv_se)),'b',1:length(deltan),deltan,'r'); grid
    figure(2); plot(1:length(deltan),1./sqrt(diag(inv_se)),'b',1:length(deltan),-1./sqrt(diag(inv_se)),'b'); grid
    figure(3); plot(dx2,'o-'); title('dx2'); grid
    figure(4); plot(dx1*dx2,'o-'); title('DX1 * dx2'); grid

    figure(5); plot(dx1*dx2,'o-'); title('DX1 * dx2'); axis([0 80 -200 +200]); grid

%    dx1 = r + k' * inv_se * k; dx1 = inv(dx1);
    figure(5); pcolor(log10(abs(dx1_0)));   colorbar; colormap jet; caxis([-5 +5]); shading flat
    figure(6); pcolor(log10(abs(dx1)));   colorbar; colormap jet; caxis([-5 +5]); shading flat
    figure(7); plot(r); colorbar; colormap jet; shading flat
    figure(7); pcolor(log10(abs(dx1 * dx1_0))); colorbar; colormap jet; shading flat
    plot(dx1)
    figure(8); plot(k); colorbar; colormap jet; shading flat
  
    pause(0.1)
    %error('lks')
  end

  % Update first guess with deltax changes
  xnbefore = xn;

  rodgers_rate = real(xn + deltax);
  figure(9); plot(1:length(xn),xn,'ko-',1:length(xn),deltax,'bx-',1:length(xn),real(xn+deltax),'r.-')
  figure(9); plot(1:length(xn),xn.*driver.qrenorm','ko-',1:length(xn),deltax.*driver.qrenorm','bx-',1:length(xn),real(xn+deltax).*driver.qrenorm','r.-')
    plotaxis2;   xlim([0 max(driver.jacobian.scalar_i)+1])
  hl = legend('orig xn','delta xn','new xn = (orig+delta)','location','best','fontsize',8); 

  ah0 = xn.*driver.qrenorm';
  dah = deltax.*driver.qrenorm';
  ah1 = real(xn+deltax).*driver.qrenorm';
  figure(8); 
    subplot(131); plot(ah0(driver.jacobian.water_i),1:length(driver.jacobian.water_i),'b',ah1(driver.jacobian.water_i),1:length(driver.jacobian.water_i),'r'); title('WV'); set(gca,'ydir','reverse')
    subplot(132); plot(ah0(driver.jacobian.temp_i),1:length(driver.jacobian.water_i),'b',ah1(driver.jacobian.temp_i),1:length(driver.jacobian.water_i),'r');  title('T'); set(gca,'ydir','reverse')
    subplot(133); plot(ah0(driver.jacobian.ozone_i),1:length(driver.jacobian.water_i),'b',ah1(driver.jacobian.ozone_i),1:length(driver.jacobian.water_i),'r'); title('O3'); set(gca,'ydir','reverse')

  xn = rodgers_rate;
  xsave(ii,:) = rodgers_rate;

  if ii <= driver.oem.nloop
    %% so this will be executed even if driver.oem.nloop == 1

    deltanIN = deltan;
    xn = rodgers_rate;

    % Form the computed rates; also see lines 108-121 of oem_lls.m
    thefitr      = zeros(1,length(driver.rateset.rates));
    thefitrdelta = zeros(1,length(driver.rateset.rates));
    for ix = 1 : length(xn)
      thefitr      = thefitr + xn(ix)*m_ts_jac(:,ix)';
      thefitrdelta = thefitrdelta + deltax(ix)*m_ts_jac(:,ix)';
    end
    figure(7); plot(f(inds),deltan,'c',f(inds),thefitrdelta(inds),'r','linewidth',2); plotaxis2; title('(c) deltaN to be fitted (r) fit')
    [~,numlay] = size(k);
    numlay = (numlay-6)/3;

    figure(8); clf; plot(f(inds),deltan,'k.-',f(inds),k(:,1:6),'linewidth',2);
      hl = legend('rate','CO2','N2O','CH4','CFC11','CFC12','stemp','location','best','fontsize',10);

    figure(9); clf; plot(f(inds),deltan,'k.-',...
                         f(inds),sum(k(:,6+0*numlay+(1:numlay)),2),f(inds),sum(k(:,6+1*numlay+(1:numlay)),2),f(inds),sum(k(:,6+2*numlay+(1:numlay)),2),'linewidth',2)
      hl = legend('rate','colWV','colT','colO3','location','best','fontsize',10);

    figure(10); clf; plot(f(inds),deltan,'k.-',f(inds),k(:,1:6),...
                         f(inds),sum(k(:,6+0*numlay+(1:numlay)),2),f(inds),sum(k(:,6+1*numlay+(1:numlay)),2),f(inds),sum(k(:,6+2*numlay+(1:numlay)),2),'linewidth',2)
      hl = legend('rate','CO2','N2O','CH4','CFC11','CFC12','stemp','colWV','colT','colO3','location','best','fontsize',10);

% keyboard_nowindow
    grid;

    iDebugNLOOP = +1;
    iDebugNLOOP = -1;
    if iDebugNLOOP > 0
      hdffile = '/home/sergio/MATLABCODE/airs_l1c_srf_tables_lls_20181205.hdf';   % what he gave in Dec 2018
      vchan2834 = hdfread(hdffile,'freq');
      f = vchan2834;
      load sarta_chans_for_l1c.mat
      f = f(ichan);

      figure(1); plot(f(inds),driver.rateset.rates(inds),'k.-',f(inds),deltan0,'kx-',f(inds),deltan,'b.-',f(inds),thefitr(inds),'g',...
		    f(inds),driver.rateset.rates(inds)-thefitr(inds)','r'); grid; title(['Rates Loop ' num2str(ii) ]);
        hl = legend('orig data','orig residual Y-f(x0)','residual at start of Nth iteration Y-f(xn-1)','f(xn)','residual still to fit','location','best');
      figure(2); plot(f(inds),driver.rateset.rates(inds)-thefitr(inds)','b',f(inds),sqrt(diag(se)),'k',f(inds),-sqrt(diag(se)),'k'); 
                 grid; title(['fit this d(spectra) next after Loop ' num2str(ii) ]);
      figure(3); plot(1:length(xb),xnIN,'b.-',1:length(xb),xn,'ro-',1:length(xb),xb,'kx-'); grid; title(['ParamsN Loop ' num2str(ii) ]);
      figure(4); plot(1:length(xb),deltax,'bo-'); grid; title(['deltax Loop ' num2str(ii) ]); 

      figure(5); clf; plot(f(inds),driver.rateset.rates(inds),'b',f(inds),deltan0,'r','linewidth',2); grid
      legend('input data Y','deltan0 = Y-f(x0)','location','best');

       figure(7); plot(f(inds),deltan); grid; figure(8); plot(dx2); grid; sum(xn-xb)

       [xnbefore(1:5) deltax(1:5) xn(1:5)]
      disp('ret'); pause
    end

    % Compute chisqr, and new deltan
    deltan  = driver.rateset.rates - thefitr';                             %% till  Jan 2021
    ddeltan = (driver.rateset.rates - tracegas_offset00) - thefitrdelta';  %% after Feb 2021
    deltan  = deltan(inds);
    ddeltan = ddeltan(inds);
    chisqr(ii) = nansum(ddeltan'.*ddeltan');
   
    iYesPlot = 1;    
    if driver.oem.doplots > 0 | iYesPlot > 0
      figure(10); clf
      figure(10); plot(f(inds),driver.rateset.rates(inds),'b.-',f(inds),driver.rateset.rates(inds) - tracegas_offset00(inds),'c.-',f(inds),thefitrdelta(inds),'k.-',f(inds),ddeltan,'r.-','linewidth',2); plotaxis2;
        hl = legend('input rates','signal''= to fit after subtracting trace gas jacs','fit','signal''-fit','location','best','fontsize',8);
        title(['rodgers.m : ADJ SIGNAL \newline obs - fit at iteration ' num2str(ii)]); pause(0.1)

      figure(11); clf
      figure(11); plot(f(inds),driver.rateset.rates(inds),'b.-',f(inds),thefitr(inds),'k.-',f(inds),deltan,'r.-','linewidth',2); plotaxis2;
        hl = legend('input rates','fit','signal-fit','location','best','fontsize',8);
        title(['rodgers.m : RAW SIGNAL \newline obs - fit at iteration ' num2str(ii)]); pause(0.1)
      %plot(f(inds),deltanIN,f(inds),deltan,'r'); 
      %title(['obs - fit at iteration ' num2str(ii)]); pause(0.1)
    end
    xnIN = xn;
  end
end

figure(10); clf
figure(10); plot(f(inds),driver.rateset.rates(inds),'b.-',f(inds),driver.rateset.rates(inds) - tracegas_offset00(inds),'c.-',f(inds),thefitrdelta(inds),'k.-',f(inds),ddeltan,'r.-','linewidth',2); plotaxis2;
        hl = legend('input rates','signal''= to fit after subtracting trace gas jacs','fit','signal''-fit','location','best','fontsize',8);
      title(['rodgers.m : ADJ SIGNAL \newline obs - fit at iteration ' num2str(ii)]); pause(0.1)

figure(11); clf
figure(11); plot(f(inds),driver.rateset.rates(inds),'b.-',f(inds),thefitr(inds),'k.-',f(inds),deltan,'r.-','linewidth',2); plotaxis2;
        hl = legend('input rates','fit','signal-fit','location','best','fontsize',8);
        title(['rodgers.m : RAW SIGNAL \newline obs - fit at iteration ' num2str(ii)]); pause(0.1)

%disp('rodgers.m 2'); keyboard_nowindow

bestloop = find(chisqr ==  min(chisqr),1);
fprintf(1,'bestloop (lowest chisqr) occured at iteration %3i \n',bestloop)
rodgers_rate = xsave(bestloop,:);

if driver.oem.nloop >= 0
  fprintf(1,'printing out successive chisqr values (upto N-1 th iterate) ... %8.6f %8.6f \n',[chisqr0 chisqr(end)])
end

iDebug = 0;   %% minimum debug
iDebug = +1;  %% tons of debug
iDebug = -1;  %% NO debug
if iDebug > 0
  %% gory detail
  renormalize = driver.qrenorm';
  wah = rodgers_rate.*renormalize';
  [xb(1:5)'; zeros(1,5);  xsave(:,1:5);  zeros(1,5);  wah(1:5)]
  disp('ret 2  to continue'); pause;
elseif iDebug == 0
  %% just the basics
  renormalize = driver.qrenorm';
  wah = rodgers_rate.*renormalize';
  wah(1:5)
end

%if iAddXB > 0   %% old wierd way of doing things
%  %%% actually we never did this
%  %%% rodgers_rate = rodgers_rate + xb;  %% if you started out with non zero z priori
%end

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
elseif invtype == 4
  errorx = inverse_ridge_regression_matrix(k' * inv_se * k + r,kmax);     %% decided pinv is too unstable, Aug 2018
elseif invtype == 5
  errorx = inverse_minimum_eigenvalue_matrix_optim(k' * inv_se * k + r,kmaxrange,sigminrange,'errorx');     %% decided pinv is too unstable, Aug 2018
elseif invtype == 9999
  AKstuff = (k' * inv_se * k + r);
  [L,U] = lu(AKstuff);
  errorx = inv(U)*inv(L);                %% trying LU
end

%{
keyboard_nowindow; 
figure(1); plot(rc); title('rc = tikonov regularization');
figure(2); plot(rcov); title('rcov = cov regularization');
figure(3); plot(k' * inv_se * k); title('k'' x se-1 x k');
%}

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
elseif invtype == 4
  inv_r       = inverse_ridge_regression_matrix(r,kmax);
  inv_r_water = inverse_ridge_regression_matrix(r_water,kmax); 
  inv_r_temp  = inverse_ridge_regression_matrix(r_temp,kmax); 
elseif invtype == 5
  inv_r       = inverse_minimum_eigenvalue_matrix_optim(r,kmaxrange,sigminrange,'inv_r');
  inv_r_water = inverse_minimum_eigenvalue_matrix_optim(r_water,kmaxrange,sigminrange,'inv_r_water'); 
  inv_r_temp  = inverse_minimum_eigenvalue_matrix_optim(r_temp,kmaxrange,sigminrange,'inv_r_temp'); 
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
  junk   = double(k * inv_r * k' + se);
    gain = inv_r *k' * inverse(junk);
  junk = double(k_water*inv_r_water*k_water'+se);
    gain_water = inv_r_water*k_water'*inverse(junk);
  junk = double(k_temp*inv_r_temp*k_temp'+se);
    gain_temp  = inv_r_temp*k_temp'*inverse(junk);
elseif invtype == 4
  gain       = inv_r *k' * inverse_ridge_regression_matrix(k * inv_r * k' + se,kmax);
  gain_water = inv_r_water*k_water'*inverse_ridge_regression_matrix(k_water*inv_r_water*k_water'+se,kmax); 
  gain_temp  = inv_r_temp*k_temp'*inverse_ridge_regression_matrix(k_temp*inv_r_temp*k_temp'+se,kmax); 
elseif invtype == 5
  gain       = inv_r *k' * inverse_minimum_eigenvalue_matrix_optim(k * inv_r * k' + se,kmaxrange,sigminrange,'gain');
  gain_water = inv_r_water*k_water'*inverse_minimum_eigenvalue_matrix_optim(k_water*inv_r_water*k_water'+se,kmaxrange,sigminrange,'gain_water'); 
  gain_temp  = inv_r_temp*k_temp'*inverse_minimum_eigenvalue_matrix_optim(k_temp*inv_r_temp*k_temp'+se,kmaxrange,sigminrange,'gain_temp'); 
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
    gain_ozone = inv_r_ozone*k_ozone'*inv(k_ozone*inv_r_ozone*k_ozone'+se);
  elseif invtype == 1 
    inv_r_ozone = pinv(r_ozone);
    gain_ozone = inv_r_ozone*k_ozone'*pinv(k_ozone*inv_r_ozone*k_ozone'+se);
  elseif invtype == 2
    inv_r_ozone = invillco(r_ozone);
    gain_ozone = inv_r_ozone*k_ozone'*invillco(k_ozone*inv_r_ozone*k_ozone'+se);
  elseif invtype == 3
    inv_r_ozone = inverse(r_ozone);
    gain_ozone = inv_r_ozone*k_ozone'*inverse(k_ozone*inv_r_ozone*k_ozone'+se);
  elseif invtype == 4
    inv_r_ozone = inverse_ridge_regression_matrix(r_ozone,kmax);
    gain_ozone = inv_r_ozone*k_ozone'*inverse_ridge_regression_matrix(k_ozone*inv_r_ozone*k_ozone'+se,kmax);
  elseif invtype == 5
    inv_r_ozone = inverse_minimum_eigenvalue_matrix_optim(r_ozone,kmaxrange,sigminrange,'inv_r_ozone');
    gain_ozone = inv_r_ozone*k_ozone'*inverse_minimum_eigenvalue_matrix_optim(k_ozone*inv_r_ozone*k_ozone'+se,kmaxrange,sigminrange,'gain_ozone');
  end
  ak_ozone = gain_ozone*k_ozone;
else
  ak_ozone = zeros(size(ak_water));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ff1 = 640; ff2 = 840;
ff1 = 640; ff2 = 1640;
%figure(10); axis([ff1 ff2 -0.2 +0.2]); grid minor; 
%figure(1);  axis([ff1 ff2 -0.2 +0.2]);  grid minor;
%figure(7);  axis([ff1 ff2 -0.2 +0.2]);  grid minor; 
%figure(11); axis([ff1 ff2 -0.2 +0.2]); grid minor;
figure(10); xlim([ff1 ff2]);
figure(1);  xlim([ff1 ff2]);
figure(7);  xlim([ff1 ff2]);
figure(11); xlim([ff1 ff2]);

