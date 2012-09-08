disp('              |      LLS     OEM')
disp('--------------|------------------')

coeffs  = driver.lls.coeffs;
coeffsR = driver.oem.coeffs;

if driver.jacobian.co2
  iAX = 1;
  array = [coeffs(iAX)  coeffsR(iAX)]*2.2;
  fprintf(1,'CO2   ppmv/yr |   %8.4f %8.4f \n',array);
  end
if driver.jacobian.o3
  iAX = 2;
  array = [coeffs(iAX)  coeffsR(iAX)]*0.01;
  fprintf(1,'O3    frac/yr |   %8.4f %8.4f \n',array);
  end
if driver.jacobian.n2o
  iAX = 3;
  array = [coeffs(iAX)  coeffsR(iAX)]*1.0;
  fprintf(1,'N2O   ppbv/yr |   %8.4f %8.4f \n',array);
  end
if driver.jacobian.ch4
  iAX = 4;
  array = [coeffs(iAX)  coeffsR(iAX)]*5.0;
  fprintf(1,'CH4   ppbv/yr |   %8.4f %8.4f \n',array);
  end
if driver.jacobian.cfc11
  iAX = 5;
  array = [coeffs(iAX)  coeffsR(iAX)]*1.0;
  fprintf(1,'CFC11 pptv/yr |   %8.4f %8.4f \n',array);
  end
if driver.jacobian.stemp
  iAX = 6;
  array = [coeffs(iAX)  coeffsR(iAX)]*0.1;
  fprintf(1,'ST       K/yr |   %8.4f %8.4f \n',array);
  end

%%% see get_jacs.m
qstjacindex = [driver.jacobian.co2    ...
               driver.jacobian.o3     ...
               driver.jacobian.n2o    ...
               driver.jacobian.ch4    ...
               driver.jacobian.cfc11  ...
               driver.jacobian.stemp];

wvjacindex = zeros(1,97);
wvjacindex(driver.jacobian.wvjacindex) = 1;

tjacindex = zeros(1,97);
tjacindex(driver.jacobian.tjacindex) = 1;

alljacindex = logical([qstjacindex(1:5) wvjacindex tjacindex qstjacindex(6)]);
%%% see get_jacs.m

lengases = sum(qstjacindex);
lenwv = sum(wvjacindex);
lenT  = sum(tjacindex);

WV_retr_ind = (1:lenwv)+ lengases;
figure(2); clf
  plot(coeffsR(WV_retr_ind)*0.01,driver.jacobian.wvjacindex,...
       'linewidth',2);  title('frac W/yr');
hl = legend('OEM','location','northwest');
set(hl,'Fontsize',10)
set(gca,'ydir','reverse')

T_retr_ind = (1:lenT)+ lengases + lenwv;
figure(3); clf
  plot(coeffsR(T_retr_ind)*0.1,driver.jacobian.tjacindex,...
       'linewidth',2);  title('frac K/yr');
hl = legend('OEM','location','northwest');
set(hl,'Fontsize',10)
set(gca,'ydir','reverse')

thedata = driver.rateset.thedata;

figure(4); clf
  plot(driver.f(inds),thedata(inds)-driver.lls.fit(inds),...
       driver.f(inds),thedata(inds)-driver.oem.fit(inds)')
hl = legend('LLS','OEM','location','northwest');
set(hl,'Fontsize',10)

figure(5); clf
  vertindex = (driver.jacobian.wvjacindex);
  subplot(121);  plot(xb(WV_retr_ind)*0.01,vertindex,'b*',...
                      coeffsR(WV_retr_ind)*0.01,vertindex,'r',...
                      erarates(WV_retr_ind)*0.1,vertindex,'k',...
                      'linewidth',2);
                 title('OEM Water frac/yr'); grid
                 hl=legend('apriori','OEM','ERA','location','west'); 
                 set(hl,'Fontsize',10);
                 set(gca,'ydir','reverse')
  vertindex = (driver.jacobian.tjacindex);
  subplot(122);  plot(xb(T_retr_ind)*0.1,vertindex,'b*',...
                      coeffsR(T_retr_ind)*0.1,vertindex,'r',...
                      erarates(T_retr_ind)*0.1,vertindex,'k',...
                      'linewidth',2);
                 hl=legend('apriori','OEM','ERA'); set(hl,'Fontsize',10);
                 title('OEM Temp K/yr'); grid
                 set(gca,'ydir','reverse')

figure(6);
  wonkb = diag(r0);     wonkb = sqrt(wonkb);
  wonkr = diag(errorx); wonkr = sqrt(wonkr);
  wonke = erastd;
  subplot(121);  
    vertindex = (driver.jacobian.wvjacindex);
    plot(wonkb(WV_retr_ind)*0.01,vertindex,'b*',...
                      wonkr(WV_retr_ind)*0.01,vertindex,'r',...  
                      wonke(WV_retr_ind)*0.01,vertindex,'k',...
                      'linewidth',2);
                 title('WV unc frac/yr'); grid
                 hl=legend('apriori','OEM','ERA','location','west'); 
                 set(hl,'Fontsize',10);
                 set(gca,'ydir','reverse')
  subplot(122);  
  vertindex = (driver.jacobian.tjacindex);
    plot(wonkb(T_retr_ind)*0.1,vertindex,'b*',...
                      wonkr(T_retr_ind)*0.1,vertindex,'r',...
                      wonke(T_retr_ind)*0.1,vertindex,'k',...
                      'linewidth',2);
                 hl=legend('apriori','OEM','ERA'); set(hl,'Fontsize',10);
                 title('T unc frac K/yr'); grid
                 set(gca,'ydir','reverse')

