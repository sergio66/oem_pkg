 %% these are common to iaSequential = -1 (one gulp) or eg [150 60 100 -1] *sequential)

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
  deltan00 = driver.rateset.rates - tracegas_offset00;    %%% << this is what we are fitting, all 2645 chans >>
  deltan   = deltan00(inds);                              %%% << this is what we are fitting, strow selected ~500 chans >>
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
nlays = length(driver.jacobian.water_i);
nTG   = length(driver.jacobian.scalar_i);
%if nnn == 66
  %%% this is 20 layers = 6 + 20 WV + 20 T + 20 Oz
  wahCO2_ST = m_ts_jac(inds,[1 nTG]);
  wahWV = m_ts_jac(inds,(1:nlays)+nTG+0*nlays);
  wahT  = m_ts_jac(inds,(1:nlays)+nTG+1*nlays);
  wahO3 = m_ts_jac(inds,(1:nlays)+nTG+2*nlays);
  figure(1); plot(f(inds),driver.rateset.rates(inds) - tracegas_offset00(inds),'kx-',f(inds),sum(wahWV'),f(inds),sum(wahT'),f(inds),sum(wahO3'),f(inds),wahCO2_ST(:,1),f(inds),wahCO2_ST(:,2),'linewidth',2); 
    plotaxis2; title('initializations :  oem\_pkg/rodgers.m : nyuk'); 
    hl = legend('signal''= to fit b-g','WVjac','Tjac','O3jac','CO2jac','STjac','location','best','fontsize',10);
  figure(1); plot(f(inds),sum(wahWV'),f(inds),sum(wahT'),f(inds),sum(wahO3'),f(inds),wahCO2_ST(:,1),f(inds),wahCO2_ST(:,2),'linewidth',2); 
    plotaxis2; title('initializations : oem\_pkg/rodgers.m : nyuk'); 
    hl = legend('WVjac','Tjac','O3jac','CO2jac','STjac','location','best','fontsize',10);
%end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    iWV = (01:nlays)+nTG; %% anomaly tile spectra
    figure(5); plot(f(inds),sum(m_ts_jac(inds,iWV),2),'.-',miaow.fout,sum(miaow.jout')*0.01); title('Should be WV jac'); xlim([640 1340]); grid; grid minor
      hl = legend('input jac','from STD/g1\_jac.mat','location','best','fontsize',10);

    if iTRPorSTD == 49  
      miaow   = load('/asl/s1/sergio/AIRSPRODUCTS_JACOBIANS/STD/temp_jac.mat');
    elseif iTRPorSTD == 1  
      miaow   = load('/asl/s1/sergio/AIRSPRODUCTS_JACOBIANS/TRP/temp_jac_new.mat');
    end
    iTz = nlays+iWV; %% anomaly tile spectra
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
