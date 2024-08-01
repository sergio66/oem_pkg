 %% these are common to iaSequential = -1 (one gulp) or eg [150 60 100 -1] *sequential)

qrenorm = driver.qrenorm;

% Jacobians
m_ts_jac = aux.m_ts_jac;
if driver.topts.dataset == 30
  %% AMSU
  junkind = 6;                       junk2ind = 1;                                   junk(:,junk2ind) = aux.m_ts_jac(:,junkind); amsu_scalar_i = junk2ind; %% stemp 
  junkind = driver.jacobian.water_i; junk2ind = junk2ind(end) + (1:length(junkind)); junk(:,junk2ind) = aux.m_ts_jac(:,junkind); amsu_water_i  = junk2ind; %% WV
  junkind = driver.jacobian.temp_i;  junk2ind = junk2ind(end) + (1:length(junkind)); junk(:,junk2ind) = aux.m_ts_jac(:,junkind); amsu_temp_i   = junk2ind; %% Tz
  m_ts_jac = junk;
  clear junk junkind junk2ind

  junk = ones(1,length(m_ts_jac));
  junkind = 6;                       junk(amsu_scalar_i) = qrenorm(junkind);
  junkind = driver.jacobian.water_i; junk(amsu_water_i)  = qrenorm(junkind);
  junkind = driver.jacobian.temp_i;  junk(amsu_temp_i)   = qrenorm(junkind);
  qrenorm = junk;
  clear junk
end

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
if driver.topts.dataset == 30
  %% AMSU
  junkind = 6;                       junk2ind = amsu_scalar_i;  junk(junk2ind) = xb(6); %% stemp
  junkind = driver.jacobian.water_i; junk2ind = amsu_water_i;   junk(junk2ind) = xb(junkind); %% WV
  junkind = driver.jacobian.temp_i;  junk2ind = amsu_temp_i;    junk(junk2ind) = xb(junkind); %% Tz
  xb = junk;
  xb = reshape(xb,length(xb),1);
  clear junk junkind junk2ind
end

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

iCommonBad = -1;
bad = find(isnan(e0));
if length(bad) > 0
  fprintf(1,'common_rodgers_initializations1.m : found %4i NaN in e0, setting to 0 \n',length(bad));
  if length(bad) >= length(e0)-20
    iCommonBad = +1;
  end
  e0(bad) = 0.0;
end
bad = find(isnan(fme));
if length(bad) > 0
  fprintf(1,'common_rodgers_initializations1.m : found %4i NaN in fme, setting to 0 \n',length(bad));
  if length(bad) >= length(fme)-20
    iCommonBad = +1;
  end
  fme(bad) = 0.0;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' >>> these are the indices where xb is not 0 ie you have initialized them <<<<')
chch = find(abs(xb) > eps); 
if length(chch) > 0
  disp(' >>> these are the indices where xb is not 0 ie you have initialized them <<<<')
  printarray([chch xb(chch) xb(chch).*qrenorm(chch)'])
end

iAddXB = -1; %% new, does this really makes more sense see eg anomaly_0dayavg_resultsXloop3try2?????
iAddXB = +1; %% orig, wierd but I think it is ok as you need raBTdelta0 = obs - tracegas_offset = obs-f(x0) = obs - f(xb)
if iAddXB > 0 & driver.topts.dataset < 30
  nyuk = find(abs(xn) > eps);
  %[nyuk xn(nyuk)]
  % Form y - F(xa), this is orig code but a little wierd!!!!!!
  tracegas_offset = zeros(size(driver.rateset.rates));
  for iy = 1 : length(xn)
     tracegas_offset = tracegas_offset + (xn(iy)*m_ts_jac(:,iy));
     if iy == length(driver.jacobian.scalar_i)
       tracegas_offset6 = tracegas_offset;
     end
     if iy == length(driver.jacobian.scalar_i) + length(driver.jacobian.water_i)
       tracegas_offsetWV = tracegas_offset - tracegas_offset6;
     end
     if iy == length(driver.jacobian.scalar_i) + length(driver.jacobian.water_i) + length(driver.jacobian.temp_i)
       tracegas_offsetT = tracegas_offset - tracegas_offset6 - tracegas_offsetWV;
     end
  end
  iJUNK = [driver.jacobian.scalar_i  driver.jacobian.water_i([1 end]) driver.jacobian.temp_i([1 end]) driver.jacobian.ozone_i([1 end])];
  disp('   scalar/WV/T/O3 xb');
  disp('   xb    |qrenorm   |  xb.*qrenorm')
  disp('------------------------------------')
  for iii = 1 : length(iJUNK)
    fprintf(1,'%8.4f | %8.4f | %8.4f \n', [xn(iJUNK(iii))   driver.qrenorm(iJUNK(iii))'  xn(iJUNK(iii)).*driver.qrenorm(iJUNK(iii))']')
    if iii == length(driver.jacobian.scalar_i) | iii == length(driver.jacobian.scalar_i)+2 | iii == length(driver.jacobian.scalar_i)+4
      disp('------------------------------------')
    end
  end
  disp('------------------------------------')

  tracegas_offset00 = tracegas_offset;
  raBTdeltan00 = driver.rateset.rates - tracegas_offset00;    %%% << this is what we are fitting, all 2645 chans >>
  raBTdeltan   = raBTdeltan00(inds);                              %%% << this is what we are fitting, strow selected ~500 chans >>
  raBTdeltan0  = raBTdeltan;

elseif iAddXB < 0 & driver.topts.dataset < 30
  tracegas_offset = zeros(size(driver.rateset.rates));
  tracegas_offset00 = tracegas_offset;
  tracegas_offset6 = tracegas_offset;
  raBTdeltan00 = driver.rateset.rates - tracegas_offset00;    %%% << this is what we are fitting, all chans >>
  raBTdeltan   = raBTdeltan00(inds);                 %%% << this is what we are fitting, selected chans >>
  raBTdeltan0  = raBTdeltan;

elseif driver.topts.dataset == 30
  tracegas_offset = zeros(size(driver.rateset.rates));
  tracegas_offset00 = tracegas_offset;
  tracegas_offset6 = tracegas_offset;
  raBTdeltan00 = driver.rateset.rates;               %%% << this is what we are fitting, all chans >>
  raBTdeltan   = raBTdeltan00(inds);                 %%% << this is what we are fitting, selected chans >>
  raBTdeltan0  = raBTdeltan;
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
elseif length(driver.rateset.rates) == 13
  f = aux.f;
else
  error('oooorrr is this AIRS 2378 or 2465 or Cris 1305 or AMSU 13?')
end

figure(12); plot(f,tracegas_offset,'b.-',f,tracegas_offset6,'r',f,m_ts_jac(:,1:3)); hl = legend('tracegas+T/WV/O3 offset','tracegas ONLY offset','CO2 jac','N2O jac','CH4 jac','location','best','fontsize',10);

figure(1); plot(f(inds),driver.rateset.rates(inds),'b.-',f(inds),tracegas_offset(inds),'g.-',f(inds),driver.rateset.rates(inds) - tracegas_offset00(inds),'c.-','linewidth',2); 
  plotaxis2; title('in oem\_pkg/rodgers.m : nyuk'); 
  hl = legend('input rates','trace gas jacs offset','signal''= to fit b-g','location','best');

if driver.topts.dataset < 30
  [mmm,nnn] = size(m_ts_jac);
  nlays = length(driver.jacobian.water_i);
  nTG   = length(driver.jacobian.scalar_i);
  %if nnn == 66
    %%% this is 20 layers = 6 + 20 WV + 20 T + 20 Oz
    wahCO2_ST = m_ts_jac(inds,[1 nTG]);
    wahWV = m_ts_jac(inds,(1:nlays)+nTG+0*nlays);
    wahT  = m_ts_jac(inds,(1:nlays)+nTG+1*nlays);
    wahO3 = m_ts_jac(inds,(1:nlays)+nTG+2*nlays);
    figure(1); plot(f(inds),driver.rateset.rates(inds) - tracegas_offset00(inds),'kx-',...
                    f(inds),sum(wahWV'),f(inds),sum(wahT'),f(inds),sum(wahO3'),f(inds),wahCO2_ST(:,1),f(inds),wahCO2_ST(:,2),'linewidth',2); 
      plotaxis2; title('initializations :  oem\_pkg/rodgers.m : nyuk'); 
      hl = legend('signal''= to fit b-g','WVjac','Tjac','O3jac','CO2jac','STjac','location','best','fontsize',10);
    figure(1); plot(f(inds),sum(wahWV'),f(inds),sum(wahT'),f(inds),sum(wahO3'),f(inds),wahCO2_ST(:,1),f(inds),wahCO2_ST(:,2),'linewidth',2); 
      plotaxis2; title('initializations : oem\_pkg/rodgers.m : nyuk'); 
      hl = legend('WVjac','Tjac','O3jac','CO2jac','STjac','location','best','fontsize',10);
    figure(1); plot(f(inds),sum(wahWV')/10,f(inds),sum(wahT'),f(inds),sum(wahO3')/10,f(inds),wahCO2_ST(:,1)/10,f(inds),wahCO2_ST(:,2),'linewidth',2); 
      plotaxis2; title('initializations : oem\_pkg/rodgers.m : nyuk'); 
      hl = legend('WVjac/10','Tjac','O3jac/10','CO2jac/10','STjac','location','best','fontsize',10);
  %end
end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if iAddXB > 0 & driver.topts.dataset < 30
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
