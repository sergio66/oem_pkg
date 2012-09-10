clear al
paper7_rates_bw
disp('ret to continue : '); pause

amsurates = input('enter (1) 2002/08-2004/08 or (2) 2002/08-2011/08  (3) 2002/08-2010/08  RSS AMSU rates : ');

if amsurates == 1
  amsufile = '/strowdata1/shared/sergio/MATLABCODE/RATES_GLOBALVIEW/amsu_rateNEW_2002_9_2004_8.mat';
elseif amsurates == 2
  amsufile = '/strowdata1/shared/sergio/MATLABCODE/RATES_GLOBALVIEW/amsu_rateNEW_2002_9_2011_8.mat';
elseif amsurates == 3
  amsufile = '/strowdata1/shared/sergio/MATLABCODE/RATES_GLOBALVIEW/amsu_rateNEW_2002_9_2010_8.mat';
else
  error('ooops amsu');
end
load(amsufile)

%% close all
load /home/sergio/MATLABCODE/RATES_TARO/MAT/overocean_gsx_1dayV1_era3_lays_spanday01_save_lat_Sep7_2011.mat

plevs = load('/home/sergio/MATLABCODE/airslevels.dat');
playsA = plevs(1:100)-plevs(2:101);
playsB = log(plevs(1:100)./plevs(2:101));
plays = playsA./playsB;
plays = plays(4:100);
plays = flipud(plays);

tropics = find(abs(save_lat) <= 30);
nml     = find(save_lat > +30  & save_lat <= 60);
sml     = find(save_lat >= -60 & save_lat < -30);
np      = find(save_lat > +60);
sp      = find(save_lat < -60);
ml      = [sml nml];
pl      = [sp  np];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% these are ERA rates
%  ozonerate          36x97            13968  single              
%  ozoneratestd       36x97            27936  double              
%  ptemprate          36x97            13968  single              
%  ptempratestd       36x97            27936  double              
%  stemprate           1x36              144  single              
%  stempratestd        1x36              288  double              
%  timeused          159x1              1272  double              
%  waterrate          36x97            13968  single              
%  waterratestd       36x97            27936  double              
% these are ERA rates

% these are OEM rates
% wv,dwv  and t,dt
% these are OEM rates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(6); clf
plot(mean(double(waterrate(tropics,:))),plays,'--','color','k','linewidth',2); 
  hold on
plot(mean(double(waterrate(nml,:))),plays,'--','color',[0.7 0.7 0.7],'linewidth',2); 
plot(mean(double(waterrate(sml,:))),plays,'--','color',[0.9 0.9 0.9],'linewidth',2); 
plot(mean(double(wv(tropics,:))),plays,'color','k','linewidth',4); 
plot(mean(double(wv(nml,:))),plays,'color',[0.7 0.7 0.7],'linewidth',4); 
plot(mean(double(wv(sml,:))),plays,'color',[0.9 0.9 0.9],'linewidth',4); 
  hold off
xlabel('WV rates (fraction/yr)')
title('Solid : OEM  Dashed : ERA');
hl = legend('tropics','N mid lat','S mid lat','location','east'); 
  set(hl,'Fontsize',10);
axis([-0.02 +0.02 100 1000]); 
set(gca,'ydir','reverse'); grid

tji = jacobian.tjacindex;
figure(7); clf
plot(10*mean(double(ptemprate(tropics,:))),plays,'--','color','k','linewidth',2); 
  hold on
plot(10*mean(double(ptemprate(nml,:))),plays,'--','color',[0.7 0.7 0.7],'linewidth',2); 
plot(10*mean(double(ptemprate(sml,:))),plays,'--','color',[0.9 0.9 0.9],'linewidth',2); 
plot(10*mean(double(t(tropics,:))),plays(tji),'color','k','linewidth',4); 
plot(10*mean(double(t(nml,:))),plays(tji),'color',[0.7 0.7 0.7],'linewidth',4); 
plot(10*mean(double(t(sml,:))),plays(tji),'color',[0.9 0.9 0.9],'linewidth',4); 
  hold off
xlabel('Temperature rates (K/decade)')
title('Solid : OEM  Dashed : ERA');
hl = legend('tropics','N mid lat','S mid lat','location','east'); 
set(hl,'Fontsize',10);
axis([-1 +1 100 1000]); 
set(gca,'ydir','reverse'); grid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(8)
ind1 = 1 : 2 : length(plays); ind1 = unique([1 ind1 length(plays)]); ix1 = 2; 
ind2 = 2 : 2 : length(plays); ind2 = unique([2 ind2 length(plays)]); ix2 = 3; 

%% ERA followed by OEM
hl=errorbar_x(mean(double(waterrate(tropics,ind1))),plays(ind1),...
              mean(waterratestd(tropics,ind1)),'o-');
  set(hl,'color',[0.5 0.5 0.5]); set(hl,'linewidth',ix1); 
hold on
hl=errorbar_x(mean(double(wv(tropics,ind2))),plays(ind2),mean(dwv(tropics,ind2)),'o-');
  set(hl,'color','k'); set(hl,'linewidth',ix2); 

hl=errorbar_x(mean(double(waterrate(nml,ind1)) + 0.02),plays(ind1),...
              mean(waterratestd(nml,ind1)),'o-');
  set(hl,'color',[0.5 0.5 0.5]); set(hl,'linewidth',ix1); 
hl=errorbar_x(mean(double(wv(nml,ind2))) + 0.02,plays(ind2),mean(dwv(nml,ind2)),'o-');
  set(hl,'color','k'); set(hl,'linewidth',ix2); 

hl=errorbar_x(mean(double(waterrate(sml,ind1))- 0.02),plays(ind1),...
              mean(waterratestd(sml,ind1)),'o-');
  set(hl,'color',[0.5 0.5 0.5]); set(hl,'linewidth',ix1); 
hl=errorbar_x(mean(double(wv(sml,ind2))-0.02),plays(ind2),mean(dwv(sml,ind2)),'o-');
  set(hl,'color','k'); set(hl,'linewidth',ix2); 
hold off
set(gca,'ydir','reverse'); grid; axis([-0.06 +0.06 100 1000])

text(-0.01,800,'TRP','fontsize',12);
text(+0.04,800,'NML','fontsize',12);
text(-0.04,800,'SML','fontsize',12);
title('WV Rates (fraction/yr) : black = OEM  gray = ERA')

tji2 = ind2;
figure(9); clf
offset = 1.0;
%% ERA followed by OEM
hl=errorbar_x(10*mean(double(ptemprate(tropics,ind1))),plays(ind1),...
              10*mean(ptempratestd(tropics,ind1)),'o-');
  set(hl,'color',[0.5 0.5 0.5]); set(hl,'linewidth',ix1); 
hold on
hl=errorbar_x(10*mean(double(t(tropics,ind2))),plays(tji2),10*mean(dt(tropics,ind2)),'o-');
  set(hl,'color','k'); set(hl,'linewidth',ix2); 

hl=errorbar_x(10*mean(double(ptemprate(nml,ind1))) + offset,plays(ind1),...
              10*mean(ptempratestd(nml,ind1)),'o-');
  set(hl,'color',[0.5 0.5 0.5]); set(hl,'linewidth',ix1); 
hl=errorbar_x(10*mean(double(t(nml,ind2))) + offset,plays(tji2),10*mean(dt(nml,ind2)),'o-');
  set(hl,'color','k'); set(hl,'linewidth',ix2); 

hl=errorbar_x(10*mean(double(ptemprate(sml,ind1)))- offset,plays(ind1),...
              10*mean(ptempratestd(sml,ind1)),'o-');
  set(hl,'color',[0.5 0.5 0.5]); set(hl,'linewidth',ix1); 
hl=errorbar_x(10*mean(double(t(sml,ind2)))-offset,plays(tji2),10*mean(dt(sml,ind2)),'o-');
  set(hl,'color','k'); set(hl,'linewidth',ix2); 
hold off
set(gca,'ydir','reverse'); grid; axis([-2*offset +2*offset 100 1000])
title('T Rates (K/decade) : black = OEM  gray = ERA')

figure(10)
% http://www.ssmi.com/msu/msu_data_description.html#msu_weighting_functions
offset = 1.0;
%% ERA followed by OEM
hl=errorbar_x(10*mean(double(ptemprate(tropics,ind1))),plays(ind1),...
              10*mean(ptempratestd(tropics,ind1)),'o-');
  set(hl,'color',[0.5 0.5 0.5]); set(hl,'linewidth',ix1); 
hold on
hl=errorbar_x(10*mean(double(t(tropics,ind2))),plays(tji2),10*mean(dt(tropics,ind2)),'o-');
  set(hl,'color','k'); set(hl,'linewidth',ix2); 
%% amsu file is LS TS MT LT = 150 300 600 800 mb
amsu_p = [150 300 600 800];
amsu_p = h2p([17000 10000 5000 2000]);
hl=errorbar_x(10*amsu_rss_rate_robust,amsu_p,10*amsu_rss_error_robust,'s--');
  set(hl,'color',[0.75 0.75 0.75]); set(hl,'linewidth',ix2); 
hold off
set(gca,'ydir','reverse'); grid; axis([-0.5 +0.5 100 1000])
title('trop T Rates (K/decade) : black = OEM  gray = ERA darkgray = AMSU RSS')

%%%%%%%%%%%%%%%%%%%%%%%%%

ppplevs = load('/home/sergio/MATLABCODE/airslevels.dat');
ppplevs = ppplevs(5:101); ppplevs = flipud(ppplevs);

figure(3)
pcolor(save_lat,ppplevs,10*double(ptemprate')); shading flat; title('ERA K/decade');
set(gca,'ydir','reverse');
caxis([-0.2 +0.2]); colorbar; colormap(usa2)
caxis([-0.5 +0.5]); colorbar; colormap(usa2)

figure(4)
pcolor(save_lat,ppplevs,10*double(t')); shading flat; title('OEM K/decade');
set(gca,'ydir','reverse');
caxis([-0.2 +0.2]); colorbar; colormap(usa2)
caxis([-0.5 +0.5]); colorbar; colormap(usa2)

%%%%%%%%%%%%%%%%%%%%%%%%%
figure(5)
pcolor(save_lat,ppplevs,double(waterrate')); shading flat; title('ERA WV frac/yr');
set(gca,'ydir','reverse');
caxis([-0.025 +0.025]); colorbar; colormap(usa2)

figure(6)
pcolor(save_lat,ppplevs,double(wv')); shading flat; title('OEM WV frac/year');
set(gca,'ydir','reverse');
caxis([-0.025 +0.025]); colorbar; colormap(usa2)

rateset

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(8)
ix = 4; ix = 1;
%% ERA
hl=errorbar_x(mean(double(waterrate(tropics,:))),plays,...
              mean(waterratestd(tropics,:)),'o--');
  set(hl,'color','k'); set(hl,'linewidth',ix); 
hold on
hl=errorbar_x(mean(double(waterrate(nml,:))),plays,...
              mean(waterratestd(nml,:)),'o--');
  set(hl,'color',[0.7 0.7 0.7]); set(hl,'linewidth',ix); 
hl=errorbar_x(mean(double(waterrate(sml,:))),plays,...
              mean(waterratestd(sml,:)),'o--');
  set(hl,'color',[0.9 0.9 0.9]); set(hl,'linewidth',ix); 
%% OEM
hl=errorbar_x(mean(double(wv(tropics,:))),plays,mean(dwv(tropics,:)),'o-');
  set(hl,'color',[0 0 0]); set(hl,'linewidth',ix); 
hl=errorbar_x(mean(double(wv(nml,:))),plays,mean(dwv(nml,:)),'o-');
  set(hl,'color',[0.7 0.7 0.7]); set(hl,'linewidth',ix); 
hl=errorbar_x(mean(double(wv(sml,:))),plays,mean(dwv(sml,:)),'o-');
  set(hl,'color',[0.9 0.9 0.9]); set(hl,'linewidth',ix); 
hold off
set(gca,'ydir','reverse'); grid; axis([-0.04 +0.04 100 1000])

figure(9)
ix = 4; ix = 1;
%% ERA
hl=errorbar_x(10*mean(double(ptemprate(tropics,:))),plays,...
              mean(ptempratestd(tropics,:)),'o--');
  set(hl,'color','k'); set(hl,'linewidth',ix); 
hold on
hl=errorbar_x(10*mean(double(ptemprate(nml,:))),plays,...
              mean(ptempratestd(nml,:)),'o--');
  set(hl,'color',[0.7 0.7 0.7]); set(hl,'linewidth',ix); 
hl=errorbar_x(10*mean(double(ptemprate(sml,:))),plays,...
              mean(ptempratestd(sml,:)),'o--');
  set(hl,'color',[0.9 0.9 0.9]); set(hl,'linewidth',ix); 
%% OEM
hl=errorbar_x(10*mean(double(t(tropics,:))),plays,mean(dt(tropics,:)),'o-');
  set(hl,'color',[0 0 0]); set(hl,'linewidth',ix); 
hl=errorbar_x(10*mean(double(t(nml,:))),plays,mean(dt(nml,:)),'o-');
  set(hl,'color',[0.7 0.7 0.7]); set(hl,'linewidth',ix); 
hl=errorbar_x(10*mean(double(t(sml,:))),plays,mean(dt(sml,:)),'o-');
  set(hl,'color',[0.9 0.9 0.9]); set(hl,'linewidth',ix); 
hold off
title('K/decade')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iL = find(plays > 400 & plays < 950);
figure(1); clf
plot(mean(double(waterrate(tropics,iL))),mean(double(wv(tropics,iL))),...
     'o','color','k','linewidth',4); 
  hold on
plot(mean(double(waterrate(nml,iL))),mean(double(wv(nml,iL))),...
     's','color',[0.7 0.7 0.7],'linewidth',4); 
plot(mean(double(waterrate(sml,iL))),mean(double(wv(sml,iL))),...
     'd','color',[0.9 0.9 0.9],'linewidth',4); 
  hold off
xlabel('ERA WV rates (fraction/yr)'); ylabel('OEM WV rates (fraction/yr)')
line([-0.04 0.04],[-0.04 0.04],'color','k')
axis([-0.02 +0.04 -0.02 +0.02]); grid

figure(2); clf
plot(10*mean(double(ptemprate(tropics,iL))),10*mean(double(t(tropics,iL))),...
     'o','color','k','linewidth',4); 
  hold on
plot(10*mean(double(ptemprate(nml,iL))),10*mean(double(t(nml,iL))),...
     's','color',[0.7 0.7 0.7],'linewidth',4); 
plot(10*mean(double(ptemprate(sml,iL))),10*mean(double(t(sml,iL))),...
      'd','color',[0.9 0.9 0.9],'linewidth',4); 
  hold off
xlabel('ERA T rates (K/decade)'); ylabel('OEM T rates (K/decade)')
line([-0.15 0.15],[-0.15 0.15],'color','k')
axis([-1.0 +1.5 -0.5 +0.5]); grid

%% printfig(6,'figure3_WV_rates','png')
%% printfig(7,'figure4_T_rates','png')
% scp figure3_WV_rates* figure4_T_rates* sergio@chard.umbc.edu:///home/sergio/PAPERS/SUBMITPAPERS/PAPER7_T_WV/FIGS/.

figure(11);
hl = errorbar_x(stemp*10,save_lat,dstemp*10,'k');
set(hl,'linewidth',2)
hold on
llS = mean(save_lat(sml)); 
llT = mean(save_lat(tropics)); llN = mean(save_lat(nml));
aaS = mean(stemp(sml)); aaT = mean(stemp(tropics)); aaN = mean(stemp(nml));
plot([aaS aaT aaN],[llS llT llN],'ko','markersize',10,'linewidth',3);
hl = errorbar_x(stemprate*10,save_lat,stempratestd*10);
set(hl,'color',[0.6 0.6 0.6])
set(hl,'linewidth',2)
aaS = mean(stemprate(sml)); aaT = mean(stemprate(tropics)); 
aaN = mean(stemprate(nml));
hl = plot([aaS aaT aaN],[llS llT llN],'x','markersize',10,'linewidth',3);
set(hl,'color',[0.6 0.6 0.6])
title('Stemp K/decade'); grid
hold off; axis([-1 +1 -80 +80])
