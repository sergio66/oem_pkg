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
% DO ALL GEOPHYSICAL VARS in ONE GULP
%---------------------------------------------------------------------------

common_rodgers_initializations1

%---------------------------------------------------------------------------

get_inv_se_rcov

%---------------------------------------------------------------------------

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
    deltan = driver.rateset.rates - thefitr';                             %% till  Jan 2021
    zeltan = (driver.rateset.rates - tracegas_offset00) - thefitrdelta';  %% after Feb 2021
    deltan = deltan(inds);
    zeltan = zeltan(inds);
    chisqr(ii) = nansum(zeltan'.*zeltan');
   
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

figure(10); clf
figure(10); plot(f(inds),driver.rateset.rates(inds),'b.-',f(inds),driver.rateset.rates(inds) - tracegas_offset00(inds),'c.-',f(inds),thefitrdelta(inds),'k.-',f(inds),zeltan,'r.-','linewidth',2); plotaxis2;
        hl = legend('input rates','signal''= to fit after subtracting trace gas jacs','fit','signal''-fit','location','best','fontsize',8);
      title(['rodgers.m : ADJ SIGNAL \newline obs - fit at iteration ' num2str(ii)]); pause(0.1)

figure(11); clf
figure(11); plot(f(inds),driver.rateset.rates(inds),'b.-',f(inds),thefitr(inds),'k.-',f(inds),deltan,'r.-','linewidth',2); plotaxis2;
        hl = legend('input rates','fit','signal-fit','location','best','fontsize',8);
        title(['rodgers.m : RAW SIGNAL \newline obs - fit at iteration ' num2str(ii)]); pause(0.1)

%[xb(1:10) xn(1:10)]
%figure(11); title('One gulp, one stage'); disp('rodgers.m 2'); keyboard_nowindow

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

do_the_dof_avg_kernel

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

