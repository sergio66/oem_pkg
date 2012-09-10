iCH4 = 1;

co2 = [];  dco2 = [];
o3  = [];  do3  = [];
n2o = [];  dn2o = [];
ch4 = [];  dch4 = [];
cfc = [];  dcfc = [];
stemp = [];dstemp = [];
fitted_rates = [];
input_rates  = [];

t = [];  dt = [];
wv = []; dwv = [];

testtype = 1;   %  L1
testtype = 0;   %  L0
testtype = 2;   %  L2

load latbins.mat

testtype = input('enter 0.1.2 for testtype : ');

%ecmrates = input('enter (1) May 23,2011 or (2) Aug 18, 2011    ERA rates : ');
ecmrates = input('enter (1) 2002-2011 or (2) 2007-2011  or (3) 2007-2012  ERA rates : ');
ecmfile = '/strowdata1/shared/sergio/MATLABCODE/RATES_TARO/MAT/';
if ecmrates == 1
  co2str = '_07_2002_07_2010';
  ecmfile = [ecmfile ...
      'overocean_gsx_1dayV1_era2_lays_profilerates_May23_2011_robust.mat'];
elseif ecmrates == 2
  %ecmfile = [ecmfile ...
  %    'overocean_gsx_1dayV1_ecmwf2_lays_spanday16_profilerates_Aug18_2011_robust.mat'];
  co2str = '_07_2007_07_2010';
  ecmfile = [ecmfile ...
      'overocean_gsx_1dayV1_era_lays_spanday01_profilerates_Aug26_2011_robust.mat'];
elseif ecmrates == 3
  co2str = '_07_2007_07_2010';
  ecmfile = [ecmfile ...
    'overocean_gsx_1day_clr_era_lays_spanday01_profilerates_Aug28_2012_robust_span_07_2007_07_2012.mat'];
end
load(ecmfile);

disp('tracegas times : (1) 8yr 07/2002-07/2010 (2) 8 yr 09/2003-09/2011 (3) 9 yr 09/2002-09/2011');
disp('tracegas times : (1) 8yr 09/2002-08/2010 (2) 7 yr 09/2003-08/2010')
tracegasrates = input('Enter tracegas rates : ');
if tracegasrates == 1
  co2str = '_07_2002_07_2010';
  co2str = '_09_2002_08_2010';
elseif tracegasrates == 2
  co2str = '_09_2003_09_2011';
  co2str = '_09_2003_08_2010';
elseif tracegasrates == 3
  co2str = '_09_2003_09_2011';
end

for iibin = 1 : 36
  clear oem
  fname = ['../Output/test' int2str(testtype) '_' int2str(iibin) '.mat'];
  fname = ['../Output/testx_' int2str(iibin) '.mat'];
  load(fname)

  qstjacindex = [jacobian.qstYesOrNo];
  qstjacindex = logical(qstjacindex);

  if length(qstjacindex) == 10
    qst10
  elseif length(qstjacindex) == 6
    qst6
  elseif length(qstjacindex) == 9
    qst9
  end

  wvjacindex = (1 : length(jacobian.Q1jacindex)) + sum(qstjacindex);
  if length(qstjacindex) ~= 9
    tjacindex  = (1 : length(jacobian.tjacindex))  + sum(qstjacindex) + ...
                  length(wvjacindex);
  elseif length(qstjacindex) == 9
    hdojacindex  = (1 : length(jacobian.tjacindex))  + sum(qstjacindex) + ...
                  length(wvjacindex);
    tjacindex  = (1 : length(jacobian.tjacindex))  + sum(qstjacindex) + ...
                  2*length(wvjacindex);
    hdo(iibin,:)  = oem.finalrates(hdojacindex); 
    dhdo(iibin,:) = oem.finalsigs(hdojacindex); 
  end

  wv(iibin,:)  = oem.finalrates(wvjacindex); 
  dwv(iibin,:) = oem.finalsigs(wvjacindex);
  t(iibin,:)   = oem.finalrates(tjacindex); 
  dt(iibin,:)  = oem.finalsigs(tjacindex);

  fitted_rates(iibin,:) = oem.fit;
  input_rates(iibin,:)  = rateset.rates;

  ak_wv = oem.ak(wvjacindex,wvjacindex);
  ak_t  = oem.ak(tjacindex,tjacindex);
  
  wv_era_smoothed(iibin,:)=ak_wv'*double(waterrate(iibin,jacobian.Q1jacindex)');
  t_era_smoothed(iibin,:) =ak_t' *double(ptemprate(iibin,jacobian.tjacindex)');
  end


if findstr(rateset.datafile,'/')
  load(rateset.datafile);   %% absolute path
else
  load(['../Test/' rateset.datafile]);  %% local path
end

switch rateset.ocb_set
  case 'bias'
     rates = squeeze(b_bias(1:36,:,2));
     unc_rates = squeeze(b_err_bias(1:36,:,2));
  case 'cal'
     rates = squeeze(b_cal(1:36,:,2));
     unc_rates = squeeze(b_err_cal(1:36,:,2));
  case {'obs','tracegas'}
     rates = squeeze(b_obs(1:36,:,2));
     unc_rates = squeeze(b_err_obs(1:36,:,2));
end

xstartup
%load ../CFC/cfcrate.mat
load ../CFC_IASI/cfcrate_IASI.mat

gray = [0.5 0.5 0.5];

figure(1); clf
hl = plot(co2,save_lat,'k');  hold on
hl = plot(n2o,save_lat);
  set(hl,'color',[0.8 0.8 0.8])
if iCH4 > 0 & jacobian.ch4
  hl = plot(ch4,save_lat);
  set(hl,'color',[0.6 0.6 0.6])
end
hl = plot(x2save36*260*0.001,save_lat);
set(hl,'color',[0.4 0.4 0.4])

if ~jacobian.co2
  disp('no CO2 in retrieval')
else
  hl = errorbar_x(co2,save_lat,dco2,'k');  hold on
  set(hl,'linewidth',2)
end

if ~jacobian.n2o
  disp('no N2O in retrieval')
else
  hl = errorbar_x(n2o,save_lat,dn2o);
  set(hl,'color',[0.8 0.8 0.8])
  set(hl,'linewidth',2); hold on
end

iCH4 = +1;
iCH4 = -1;
if ~jacobian.ch4
  disp('no CH4 in retrieval')
elseif iCH4 > 0 & jacobian.ch4
  hl = errorbar_x(ch4,save_lat,dch4);
  set(hl,'color',[0.6 0.6 0.6])
  set(hl,'linewidth',2);; hold on
end

%%errorbar_x(cfc,save_lat,dcfc,'c'); hold on
%hl = errorbar_x(x2save36*260,save_lat,save_error*260);
%set(hl,'color',[0.4 0.4 0.4])
%set(hl,'linewidth',2); 
%hold on
if ~jacobian.cfc
  disp('no CFC in retrieval')
else
  hl = errorbar_x(cfc,save_lat,dcfc);
  set(hl,'color',[0.4 0.4 0.4])
  set(hl,'linewidth',2);; hold on
end

%if ~jacobian.stemp
%  disp('no STEMP in retrieval')
%else
%  errorbar_x(stemp*10,save_lat,dstemp*10,'k'); hold on
%end
plot_globalview2_bw(co2str,iCH4);
grid

if iCH4 < 0
  %title('(B)CO2 (G)N2O (C)CFC11','fontsize',10)
  title('Rates per year','fontsize',10)
  hl = legend('CO2 (ppm)','N2O (ppb)','CFC11 (ppt)','Location','north');
  set(hl,'Fontsize',10)
  axis([-8 +3 -90 +90]);
else
  %title('(B)CO2 (G)N2O (R)CH4 (C)CFC11','fontsize',10)
  title('Rates per year','fontsize',10)
  hl = legend('CO2 (ppm)','N2O (ppb)','CH4 (ppb)','CFC11 (ppt)','Location','north');
  set(hl,'Fontsize',10)
  axis([-8 +10 -90 +90]);
end

xlabel('Rate in ppx/yr'); ylabel('latitude')
 axis([-10 +5 -90 +90]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if length(qstjacindex) > 6
  figure(8);
  plot(co2,save_lat,co2strat,save_lat,'b--',n2o,save_lat,ch4,save_lat,cfc,save_lat,'linewidth',2)
    hl = legend('CO2 (ppm)','CO2strat','N2O (ppb)','CH4 (ppb)','CFC11 (ppt)','Location','north');
  set(hl,'Fontsize',10); grid
  axis([-10 +5 -90 +90]);
end

if length(qstjacindex) == 10
  figure(9);
  plot(o3,save_lat,o3strat,save_lat,'b--',co,save_lat,hdo,save_lat,stemp,save_lat,'linewidth',2)
    hl = legend('O3','O3strat','CO','hd0','stemp','Location','north');
  set(hl,'Fontsize',10); grid
  axis([-0.03 +0.03 -90 +90]);

elseif length(qstjacindex) == 9
  figure(9);
  plot(o3,save_lat,o3strat,save_lat,'b--',co,save_lat,stemp,save_lat,'linewidth',2)
    hl = legend('O3','O3strat','CO','stemp','Location','north');
  set(hl,'Fontsize',10); grid
  axis([-0.03 +0.03 -90 +90]);

  figure(10)
  plot(wv,1:97,'b',hdo,1:97,'r'); set(gca,'ydir','reverse');
  title('(b) H2O  (r) HDO frac/yr')

  figure(11)
  plot(t,1:97,'b'); set(gca,'ydir','reverse');
  title('(b) temp Kyr')

elseif length(qstjacindex) == 6
  figure(9); clf

  figure(10)
  plot(wv,1:97,'b'); set(gca,'ydir','reverse');
  title('(b) H2O frac/yr')

  figure(11)
  plot(t,1:97,'b'); set(gca,'ydir','reverse');
  title('(b) temp Kyr')
end
%https://climatedataguide.ucar.edu/guidance/water-isotopes-satellites
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plevs = load('/home/sergio/MATLABCODE/airslevels.dat');
playsA = plevs(1:100)-plevs(2:101);
playsB = log(plevs(1:100)./plevs(2:101));
plays = playsA./playsB;
plays = plays(4:100);
plays = flipud(plays);

load /home/sergio/MATLABCODE/RATES_TARO/MAT/overocean_gsx_1dayV1_era3_lays_spanday01_save_lat_Sep7_2011.mat
tropics = find(abs(save_lat) <= 30);

figure(3);
  bwah = size(hdo);
  if bwah(1) == 97 | bwah(2) == 97
    plot(wv(18,:),plays,hdo(18,:),plays,'c',nanmean(waterrate(tropics,:)),plays,'r--','linewidth',2);
    set(gca,'ydir','reverse'); grid; title('WV frac/yr')
    hl = legend('OEM WV','OEM HDO','ERA','location','east'); set(hl,'fontsize',10)
  else
    plot(wv(18,:),plays,nanmean(waterrate(tropics,:)),plays,'r--','linewidth',2);
    set(gca,'ydir','reverse'); grid; title('WV frac/yr')
    hl = legend('OEM','ERA','location','east'); set(hl,'fontsize',10)
  end
hold on
plot(nanmean(waterrate(tropics,:))+nanstd(waterrate(tropics,:)),plays,'m--','linewidth',2);
plot(nanmean(waterrate(tropics,:))-nanstd(waterrate(tropics,:)),plays,'m--','linewidth',2);
hold off

figure(4)
plot(t(18,:),plays,nanmean(ptemprate(tropics,:)),plays,'r--','linewidth',2);
hold on
plot(nanmean(ptemprate(tropics,:))+nanstd(ptemprate(tropics,:)),plays,'m--','linewidth',2);
plot(nanmean(ptemprate(tropics,:))-nanstd(ptemprate(tropics,:)),plays,'m--','linewidth',2);
hold off
  set(gca,'ydir','reverse'); grid; title('T K/yr')
  hl = legend('OEM','ERA','location','east'); set(hl,'fontsize',10)

chanset = jacobian.chanset;
%g = dogoodchan;
figure(5);
plot(f,input_rates(18,:),'b',f,fitted_rates(18,:),'r',...
      f,input_rates(18,:)-fitted_rates(18,:),'k',...
      f(chanset),input_rates(18,chanset)-fitted_rates(18,chanset),'ko')
  hl = legend('input','fit','bias=input-fit','location','north'); set(hl,'fontsize',10)
  axis([640 2780 -0.10 +0.10]); grid

wvrates.oem = wv(18,:);
wvrates.oem_d = dwv(18,:);
wvrates.era = nanmean(waterrate(tropics,:));
wvrates.plays = plays;

trates.oem = t(18,:);
trates.oem_d = dt(18,:);
trates.era = nanmean(ptemprate(tropics,:));
trates.plays = plays;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(1,'rateset = %s \n',rateset.datafile);
fprintf(1,'jacset  = %s \n',jacobian.filename);

