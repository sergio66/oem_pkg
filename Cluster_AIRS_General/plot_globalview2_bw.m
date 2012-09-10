function boo = plot_globalview2_bw(datasavestr0,iCH4);

datasavestr = '_11_2003_11_2008';
loader = ['load ../../RATES_GLOBALVIEW/agage_rate' datasavestr '.mat'];  
  eval(loader);
loader = ...
  ['load ../../RATES_GLOBALVIEW/noaahats_cfc11_rate' datasavestr '.mat'];

datasavestr = '_07_2002_07_2010';
loader = ['load ../../RATES_GLOBALVIEW/co2_gv_rate' datasavestr '.mat']; 
  eval(loader);

datasavestr = '_11_2003_11_2008';
loader = ['load ../../RATES_GLOBALVIEW/ch4_gv_rate' datasavestr '.mat']; 
  eval(loader);

loader = ['load ../../RATES_GLOBALVIEW/co2_gv_rate' datasavestr0 '.mat']; 
  eval(loader);

bablue = [0 0 0.5];
gcf; hold on
h = plot(cfc11_agage_rate,agage_latitude,'color',[0.4 0.4 0.4],...
           'LineStyle','o',...
           'markeredgecolor',[0.4 0.4 0.4],...
           'markerfacecolor',[0.4 0.4 0.4],...
           'markersize',10);
%get(h)
hl = errorbar_x(cfc11_agage_rate,agage_latitude,cfc11_agage_error);
set(hl,'color',[0.4 0.4 0.4])
%set(hl,'linestyle','o')
hold on

plot(n2o_agage_rate,agage_latitude,'color',[0.8 0.8 0.8],...
           'LineStyle','o',...
           'markeredgecolor',[0.8 0.8 0.8],...
           'markerfacecolor',[0.8 0.8 0.8],...
           'markersize',10)
hl=errorbar_x(n2o_agage_rate,agage_latitude,n2o_agage_error);
set(hl,'color',[0.8 0.8 0.8])
%set(hl,'linestyle','o')

if iCH4 > 0
  c1 = find(ch4_gv_hgt > 5000);
  plot(ch4_gv_rate(c1),ch4_gv_latitude(c1),'color',[0.6 0.6 0.6],...
           'LineStyle','o',...
           'markeredgecolor',[0.6 0.6 0.6],...
           'markerfacecolor',[0.6 0.6 0.6],...
           'markersize',10)
  hl=errorbar_x(ch4_gv_rate(c1),ch4_gv_latitude(c1),ch4_gv_error(c1));
  %set(hl,'linestyle','o')
  set(hl,'color',[0.6 0.6 0.6])
end

c2 = find(co2_gv_hgt > 5000);
plot(co2_gv_rate(c2),co2_gv_latitude(c2),'color',[0.0 0.0 0.0],...
           'LineStyle','o',...
           'markeredgecolor',[0.0 0.0 0.0],...
           'markerfacecolor',[0.0 0.0 0.0],...
           'markersize',10)
hl = errorbar_x(co2_gv_rate(c2),co2_gv_latitude(c2),co2_gv_error(c2));
  fprintf(1,'global GV mean rate for CO2 = %8.6f +/- %8.6f ppmv/yr \n',nanmean(co2_gv_rate(c2)),nanstd(co2_gv_rate(c2)))
set(hl,'color','k')
%set(hl,'linestyle','o')

hold on; plot([-4 +6],[0 0],'k','linewidth',2); hold off;
hold off

