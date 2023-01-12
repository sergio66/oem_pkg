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
% DO sequential retrieval of GEOPHYSICAL VARS
%---------------------------------------------------------------------------

common_rodgers_initializations1

%---------------------------------------------------------------------------

iaSequential = driver.iaSequential;

seSave     = se;
rcovSave   = driver.oem.cov;
kSave      = k;
deltanSave = deltan; %% from common_rodgers_initializatio
                     %%    deltan00 = driver.rateset.rates - tracegas_offset00;    %%% << this is what we are fitting, all 2645 chans >>
                     %%    deltan   = deltan00(inds);                              %%% << this is what we are fitting, strow selected ~500 chans >>
xbSave     = xb;
xnSave     = xn;
fSave      = f;
xsave      = zeros(size(xb));

for iiS = 1 : length(iaSequential)
  iSequential = iaSequential(iiS);
  fprintf(1,' iiS = %2i of %2i : doing iSequential = %2i \n',iiS,length(iaSequential),iSequential);

  fuse = f(inds);

  if iSequential == -1
    %% nothing to do, use all chans, retrieve all geophys params
    iUseRetrParam  = 1 : length(xbSave);
    iUseChan       = 1 : length(inds);

  elseif iSequential == 150
    %% 15 um T(z) and STEMP
    iUseRetrParam = [6 driver.jacobian.temp_i];

    iUseChans15 = find(fuse <= 800);
    iUseWindow  = find((fuse > 820 & fuse <= 960) | (fuse > 1225 & fuse <= 1235));
    iUseChan    = union(iUseChans15,iUseWindow);
  elseif iSequential == 100
    %% 10 um O3(z)
    iUseRetrParam = [driver.jacobian.ozone_i];

    iUse = find(fuse > 1000 & fuse < 1080);
    iUseChan = union(iUseChans15,iUseWindow);

  elseif iSequential == 60
    %% 6 um WV(z)
    iUseRetrParam = [driver.jacobian.water_i];

    iUseWindow = find((fuse > 820 & fuse <= 960) | (fuse > 1090 & fuse <= 1280));
    iUseChans6 = find(fuse > 1380 & fuse <= 1700);
    iUseChan   = union(iUseChans6,iUseWindow);

  end

  se = seSave(iUseChan,iUseChan);
  k = kSave(iUseChan,iUseRetrParam);
  xn = xnSave(iUseRetrParam);
  xb = xbSave(iUseRetrParam);
  deltan = deltanSave(iUseChan);

  figure(21);
    plot(fuse,driver.rateset.rates(inds),'b.',fuse(iUseChan),driver.rateset.rates(inds(iUseChan)),'r')

  disp('starting iSequential loop [xbSave(1:10) xnSave(1:10) xb(1:10) xn(1:10)]')
  [xbSave(1:10) xnSave(1:10) xb(1:10) xn(1:10)]

  %-------------------------

  get_inv_se_rcov
  rcov = rcovSave(iUseRetrParam,iUseRetrParam);
  r    = r(iUseRetrParam,iUseRetrParam);

  %-------------------------
  
  % whos seSave rCovSave kSave deltanSave xbSave fuse iUseChan r k inv_se rcov rc
  chisqr0 = nansum(deltan'.*deltan');
  
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
    figure(7); plot(fuse(iUseChan),deltan); plotaxis2;                 title('deltaBT to fit')
    figure(8); plot(deltax.*driver.qrenorm(iUseRetrParam)'); plotaxis2; grid minor;   title('deltax.*qrenorm')
  
    iDebug = +1;
    iDebug = -1;
    if iDebug > 0

      %addpath /home/sergio/MATLABCODE; keyboard_nowindow
      figure(1); plot(f(inds(iUseChan)),k); grid
      figure(1); plot(f(inds(iUseChan)),k(:,1:5)); grid
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
    figure(9); plot(1:length(xn),xn.*driver.qrenorm(iUseRetrParam)','ko-',1:length(xn),deltax.*driver.qrenorm(iUseRetrParam)','bx-',1:length(xn),real(xn+deltax).*driver.qrenorm(iUseRetrParam)','r.-')
      plotaxis2;   xlim([0 max(driver.jacobian.scalar_i)+1])
    hl = legend('orig xn','delta xn','new xn = (orig+delta)','location','best','fontsize',8); 
  
    ah0 = xn.*driver.qrenorm(iUseRetrParam)';
    dah = deltax.*driver.qrenorm(iUseRetrParam)';
    ah1 = real(xn+deltax).*driver.qrenorm(iUseRetrParam)';
    figure(8); 
    if iSequential == -1
      subplot(131); plot(ah0(driver.jacobian.water_i),1:length(driver.jacobian.water_i),'b',ah1(driver.jacobian.water_i),1:length(driver.jacobian.water_i),'r'); title('WV'); set(gca,'ydir','reverse')
      subplot(132); plot(ah0(driver.jacobian.temp_i),1:length(driver.jacobian.water_i),'b',ah1(driver.jacobian.temp_i),1:length(driver.jacobian.water_i),'r');  title('T'); set(gca,'ydir','reverse')
      subplot(133); plot(ah0(driver.jacobian.ozone_i),1:length(driver.jacobian.water_i),'b',ah1(driver.jacobian.ozone_i),1:length(driver.jacobian.water_i),'r'); title('O3'); set(gca,'ydir','reverse')
    elseif iSequential == 60
      subplot(131); plot(ah0,1:length(driver.jacobian.water_i),'b',ah1,1:length(driver.jacobian.water_i),'r'); title('WV'); set(gca,'ydir','reverse')
    elseif iSequential == 150
      subplot(132); plot(ah0,1:length(ah0),'b',ah1,1:length(ah0),'r');  title('T'); set(gca,'ydir','reverse')
    elseif iSequential == 100
      subplot(133); plot(ah0,1:length(driver.jacobian.ozone_i),'b',ah1,1:length(driver.jacobian.ozone_i),'r'); title('O3'); set(gca,'ydir','reverse')
    end

    xn = rodgers_rate;
    xsave(ii,iUseRetrParam) = rodgers_rate;
  
    if ii <= driver.oem.nloop
      %% so this will be executed even if driver.oem.nloop == 1
  
      deltanIN = deltan;
      xn = rodgers_rate;
  
      % Form the computed rates; also see lines 108-121 of oem_lls.m
      thefitr      = zeros(1,length(driver.rateset.rates));
      thefitrdelta = zeros(1,length(driver.rateset.rates));

      yn = xbSave;
      yn(iUseRetrParam) = xn; %% this is needed for next iSequential

      for ix = 1 : length(xbSave)
        %% this uses all params; this uses all params;
        thefitr = thefitr + yn(ix)*m_ts_jac(:,ix)';
      end
      for ix = 1 : length(deltax)
        %% this just uses the params used in this iSequential round
        thefitrdelta = thefitrdelta + deltax(ix)*m_ts_jac(:,iUseRetrParam(ix))';
      end
      figure(7); plot(f(inds(iUseChan)),deltan,'c',f(inds),thefitrdelta(inds),'r','linewidth',2); plotaxis2; title('(c) deltaN to be fitted (r) fit')
    
      % Compute chisqr, and new deltan
      deltanIn   = deltan;
      deltan = driver.rateset.rates - thefitr';                             %% till  Jan 2021
      zeltan = (driver.rateset.rates - tracegas_offset00) - thefitrdelta';  %% after Feb 2021
      zeltan = (driver.rateset.rates - 1*tracegas_offset00) - thefitr';     %% after Jan 2023, remember we are assuming CO2/N2O/CH4 is fixed and all problems are with T/WV so should we take this out??? lets see
      zeltan = (driver.rateset.rates - 0*tracegas_offset00) - thefitr';     %% after Jan 2023, remember we are assuming CO2/N2O/CH4 is fixed and all problems are with T/WV so should we take this out??? nah
      deltanSave = deltan(inds);      %% being silly, but works beter? nah cant be
      deltanSave = zeltan(inds);      %% see eg common_rodgers_initializations1.m, note zeltan is different here than in rodgers.m
      deltan = deltan(inds);
      zeltan = zeltan(inds);
      chisqr(ii) = nansum(zeltan'.*zeltan');
     
  figure(22);
    plot(fuse,driver.rateset.rates(inds),'b-',fuse(iUseChan),driver.rateset.rates(inds(iUseChan)),'k.-',fuse(iUseChan),deltanIn,'g',fuse,deltanSave,'r')
      plotaxis2; hl = legend('all Strow chans rates','Sequential chans rates','starting delta','ending delta','location','best','fontsize',8);
  figure(23);
    if iSequential == -1
      junkST  = k(:,6);
      junkWV  = k(:,driver.jacobian.water_i);
      junkTz  = k(:,driver.jacobian.temp_i);
      junkO3  = k(:,driver.jacobian.ozone_i);
      plot(fuse(iUseChan),junkST,fuse(iUseChan),sum(junkWV,2),fuse(iUseChan),sum(junkTz,2),fuse(iUseChan),sum(junkO3,2)); plotaxis2; hl = legend('ST','WV','T','O3','location','best','fontsize',10); title('JAC')
    elseif iSequential == 150
      junkST  = k(:,1);
      junkTz  = k(:,2:length(driver.jacobian.temp_i)+1);
      plot(fuse(iUseChan),junkST,fuse(iUseChan),sum(junkTz,2)); plotaxis2; hl = legend('ST','T','location','best','fontsize',10); title('JAC')
    elseif iSequential == 100
      junkO3  = k;
      plot(fuse(iUseChan),sum(junkO3,2)); plotaxis2; hl = legend('O3','location','best','fontsize',10); title('JAC')
    elseif iSequential == 60
      junkWV  = k;
      plot(fuse(iUseChan),sum(junkWV,2)); plotaxis2; hl = legend('WV','location','best','fontsize',10); title('JAC')
    end

      iYesPlot = 1;    
      if driver.oem.doplots > 0 | iYesPlot > 0
        figure(10); clf
        figure(10); plot(f(inds),driver.rateset.rates(inds),'b.-',f(inds),driver.rateset.rates(inds) - tracegas_offset00(inds),'c.-',f(inds),thefitrdelta(inds),'k.-',f(inds),zeltan,'r.-','linewidth',2); plotaxis2;
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

  xbSave(iUseRetrParam) = xn; %% this is needed for next iSequential
  xnSave(iUseRetrParam) = xn; %% this is needed for next iSequential

  figure(10); clf
  figure(10); plot(f(inds),driver.rateset.rates(inds),'b.-',f(inds),driver.rateset.rates(inds) - tracegas_offset00(inds),'c.-',f(inds),thefitrdelta(inds),'k.-',f(inds),zeltan,'r.-','linewidth',2); plotaxis2;
          hl = legend('input rates','signal''= to fit after subtracting trace gas jacs','fit','signal''-fit','location','best','fontsize',8);
        title(['rodgers.m : ADJ SIGNAL \newline obs - fit at iteration ' num2str(ii)]); pause(0.1)
  
  figure(11); clf
  figure(11); plot(f(inds),driver.rateset.rates(inds),'b.-',f(inds),thefitr(inds),'k.-',f(inds),deltan,'r.-','linewidth',2); plotaxis2;
          hl = legend('input rates','fit','signal-fit','location','best','fontsize',8);
          title(['rodgers.m : RAW SIGNAL \newline obs - fit at iteration ' num2str(ii)]); pause(0.1)

  %disp('exiting iSequential loop [xbSave(1:10) xnSave(1:10) xb(1:10) xn(1:10)]')
  %[xbSave(1:10) xnSave(1:10) xb(1:10) xn(1:10)]
  %figure(11); title('Sequential, many stages'); disp('rodgers_sequential.m 2'); keyboard_nowindow
  
  bestloop = find(chisqr ==  min(chisqr),1);
  fprintf(1,'bestloop (lowest chisqr) occured at iteration %3i \n',bestloop)
  rodgers_rate = xsave(bestloop,:);
  
  if driver.oem.nloop >= 0
    fprintf(1,'printing out successive chisqr values (upto N-1 th iterate) ... %8.6f %8.6f \n',[chisqr0 chisqr(end)])
  end

  if iSequential == -1    
    do_the_dof_avg_kernel
  end
  
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
  
