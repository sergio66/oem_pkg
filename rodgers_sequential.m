function [rodgers_rate,errorx,dofs,cdofs,gain,ak,r,se,inv_se,se_errors,ak_water,ak_temp,ak_ozone,bestloop,raBTdeltan00] = rodgers(driver,aux)

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
%   raBTdeltan00       = ratesIN - jac*tracegas settings (spectral BT rate to be fitted)
%
% DO sequential retrieval of GEOPHYSICAL VARS
%---------------------------------------------------------------------------

common_rodgers_initializations1

%---------------------------------------------------------------------------
addpath /home/sergio/MATLABCODE/CRODGERS_FAST_CLOUD/

% function [airsL2_chans_set,airsL2_list_set,airsL2_IDs_set,fairs0,chanIDairs] = airsL2_retrieval_channel_set(nchan2378or2645,iIgnoreWindow_forWVretr,iv6or7);
% input
%   if nchan2378or2645      = +1 then we are using 2378 chan set  so eg airsL2_IDs_set = 1291,airsL2_chans_set = 1231.3
%                           = -1 then we are using 2645 chan set
%   iv6or7                  = v6 or 7 (default 7)
%   iIgnoreWindow_forWVretr = +1 (ignore when setting WV chans) or -1 (include)
% output
%    airsL2_chans_set = channel center freq    h.vchan(chosen)
%    airsL2_list_set  = channel list           chosen
%    airsL2_IDs_set   = channel ID             h.ichan(chosen)
%    fairs,chanIDairs = channels centers,cIDlist (for 2378 chans cIDlist is ordered 1:2378 ... for 2645 chans, f is sorted)
%------------------------------------------------------------------------

figure(12); clf

iaSequential = driver.iaSequential;

iaSequential_orig = iaSequential;
if iaSequential(end) ~= -1
  %% need this to do the DOF calcs 
  %% see  ~/MATLABCODE/CRODGERS_FAST_CLOUD/RODGERS/RODGERS21/rodgers21_reg_and_cov_regularization.m
  disp(' ... artificially adding on iSequential = -1 at the end, to do DofF')
  iaSequential = [iaSequential -1];
end

if (xb(1:3)) == [0 0 0] & driver.topts.ocb_set == 0
  iFitTraceGas = +1;
else
  iFitTraceGas = -1;
end
get_inv_se_rcov_allchans_allparams   %% iaSequential = -1

qrenormSave    = driver.qrenorm;
seSave         = se;
inv_seSave     = inv_se;
rcovSave       = driver.oem.cov;
rSave          = r;
%rcovSave2      = rcov2;
%rSave2         = r2;
%if exist('rcov210')
%  rcovSave210    = rcov210;
%  rSave210       = r210;
%end
%if exist('rcov214')
%  rcovSave214    = rcov214;
%  rSave214       = r214;
%end
kSave          = k;
xbSave         = xb;
xnSave         = xn;
fSave          = f;
xsave          = zeros(length(iaSequential),driver.oem.nloop,length(xb));
raBTdeltanSave = raBTdeltan; %% from common_rodgers_initialization, this is strow selected ~500 chans
                     %%    raBTdeltan00 = driver.rateset.rates - tracegas_offset00;    %%% << this is what we are fitting, all 2645 chans >>
                     %%    raBTdeltan   = raBTdeltan00(inds);                          %%% << this is what we are fitting, strow selected ~500 chans >>

raBTdeltaIterate(:,1) = raBTdeltanSave;

[airsL2_chans_set,airsL2_list_set,airsL2_IDs_set,fairs0,chanIDairs] = airsL2_retrieval_channel_set(-1,-1,7);
% keyboard_nowindow
%% PLEASE NOTE BENE that for AIRS TRENDS length(f) = 2645, and f = h.vchan for AIRS L1C so eg f(1520) = 1231.3 cm-1
%%    so basically chanIDairs(1520) = 1291 and f = fairs0

%iDebug = +1;
iDebug = -1;

iKeyBoard = +1;
iKeyBoard = -1;

fprintf(1,'iSequential = %3i : GND-2/GN-1/GND   WV and T, ST entering OEM \n',-9999);
  junk = [xbSave(driver.jacobian.water_i(end-3):driver.jacobian.water_i(end)); xbSave(driver.jacobian.temp_i(end-3):driver.jacobian.temp_i(end)); xbSave(6)];
  fprintf(1,'%8.6f %8.6f %8.6f %8.6f      %8.6f %8.6f %8.6f %8.6f     %8.6f \n',junk);

for iiS = 1 : length(iaSequential)
  disp(' ')
  iSequential = iaSequential(iiS);
  if iiS <= length(iaSequential_orig)
    iYesThisIsFine = +1;
    fprintf(1,' iiS = %2i of %2i : doing iSequential = %2i \n',iiS,length(iaSequential),iSequential);
  else
    iYesThisIsFine = -1;
    iYesThisIsFine = +1; %% well if it improves things USE IT
    fprintf(1,' iiS = %2i of %2i : need DofF so doing extra iSequential = %2i \n',iiS,length(iaSequential),iSequential);
  end

  fuse = f(inds);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  iNXYZLay = 4; %% lowest 4
  iNXYZLay = 6; %% lowest 6

  [iUseRetrParam,iUseChan] = get_sequential_indices_params(airsL2_chans_set,airsL2_list_set,airsL2_IDs_set,fairs0,chanIDairs,fuse,inds,xbSave,driver,iSequential,iNXYZLay,iFitTraceGas);
  figure(21); plot(iUseChan,f(inds(iUseChan)),'x-'); xlabel('iUseChan'); ylabel('f(inds(iUseChan))'); title(['iSequential = ' num2str(iSequential)]); 
  %disp('ret to continue'); pause

  se = seSave(iUseChan,iUseChan);          %% iSequential ~100 channels
  inv_se = inv_seSave(iUseChan,iUseChan);  %% iSequential ~100 channels
  k = kSave(iUseChan,iUseRetrParam);       %% iSequential ~100 channels,~20 params
  xn = xnSave(iUseRetrParam);              %% iSequential ~20 params
  xb = xbSave(iUseRetrParam);              %% iSequential ~20 params
  raBTdeltan = raBTdeltanSave(iUseChan);   %% iSequential ~100 channels
  qrenorm = qrenormSave(iUseRetrParam);    %% iSequential ~20 params

  get_inv_se_rcov_iSequential   %% this would make afresh : rc,rcov and thus r

  disp(' .................... >>>>>>>>>>> -------------------- <<<<<<<<<<<<< .....................')
  disp('starting iSequential loop physical vars : displaying ')
  disp('   [xbSave(ImportantParam) xnSave(ImportantParam) xb(ImportantParamX) xn(ImportantParamX)] ')
  disp(' where iImportantParam === FIXED = [CO2 N2O CH4 STEMP     WV(GND/TOA)     TZ(GND/TOA)      O3(GND/TOA) ....')
  iImportantParam  = [1 2 3 6 driver.jacobian.water_i(end) driver.jacobian.water_i(1) driver.jacobian.temp_i(end) driver.jacobian.temp_i(1) driver.jacobian.ozone_i(end) driver.jacobian.ozone_i(1)];
  if iSequential == -1
    disp(' ... and where (variable) ImportantParamX = CO2 CH4 N2O ST wvz(GND TOA) Tz(GND TOA) o3z(GND TOA)')
    iImportantParamX = [1 2 3 6    driver.jacobian.water_i(end) driver.jacobian.water_i(1) driver.jacobian.temp_i(end) driver.jacobian.temp_i(1) driver.jacobian.ozone_i(end) driver.jacobian.ozone_i(1)];
  elseif iSequential == 150
    disp(' ... and where (variable) ImportantParamX = ST ST ST ST Tz(GND TOA) Tz(GND TOA) Tz(GND TOA)')
    if iFitTraceGas == -1
      iImportantParamX = [1 1 1 1    1+length(driver.jacobian.water_i) 2 1+length(driver.jacobian.water_i) 2 1+length(driver.jacobian.water_i) 2];
    else
      iImportantParamX = [2 2 2 2    2+length(driver.jacobian.water_i) 3 2+length(driver.jacobian.water_i) 3 2+length(driver.jacobian.water_i) 3];
    end
  elseif iSequential == 210
    disp(' ... and where (variable) ImportantParamX = ST ST ST ST WVz(GND GND+4) WVz(GND GND+4) Tz(GND TOA)')
    iImportantParamX = [1 1 1 1    2+iNXYZLay-1 2 2+iNXYZLay-1 2 length(iUseRetrParam) 2+iNXYZLay];
  elseif iSequential == 214
    disp(' ... and where (variable) ImportantParamX = ST ST ST ST WVz(GND GND+4) WVz(GND GND+4) Tz(GND GND+4)')
    iImportantParamX = [1 1 1 1    2+iNXYZLay-1 2 2+iNXYZLay-1 2 length(iUseRetrParam) 2+iNXYZLay;]
  elseif iSequential == 60 | iSequential == 100
    disp(' ... and where (variable) ImportantParamX = 4 times wvz/o3z(1)     wvz/o3z(GND TOA) wvz/o3z(1) wvz/o3z(GND TOA) wvz/o3z(GND TOA) wvz/o3z(GND TOA)')
    iImportantParamX = [1 1 1 1    length(driver.jacobian.water_i) 1    length(driver.jacobian.water_i) 1    length(driver.jacobian.water_i) 1];
  end

  fprintf(1,'iSequential = %3i length(iUseRetrParam) = %3i length(iUseChan) = %3i \n',iSequential,length(iUseRetrParam),length(iUseChan))

  junk = [xbSave(iImportantParam).*driver.qrenorm(iImportantParam)' xnSave(iImportantParam).*driver.qrenorm(iImportantParam)' ...
          xb(iImportantParamX).*driver.qrenorm(iUseRetrParam(iImportantParamX))' xn(iImportantParamX).*driver.qrenorm(iUseRetrParam(iImportantParamX))'];
  for iijunk = 1 : 10
    fprintf(1,'%8.5f   %8.5f      %8.5f   %8.5f \n',junk(iijunk,:));
    if iijunk == 4 | iijunk == 6 | iijunk == 8
      disp('-------------------------------------------------------------')
    end
  end
  %-------------------------

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  figure(21);
    plot(fuse,driver.rateset.rates(inds),'b.',fuse(iUseChan),driver.rateset.rates(inds(iUseChan)),'r')

  %% get_inv_se_rcov_allchans_allparams 
  %% rcov = rcovSave(iUseRetrParam,iUseRetrParam);
  %% r    = rSave(iUseRetrParam,iUseRetrParam);
  %% if iSequential == -1
  %%   rcov = rcovSave2(iUseRetrParam,iUseRetrParam);
  %%   r    = rSave2(iUseRetrParam,iUseRetrParam);
  %% elseif iSequential == 210
  %%   rcov = rcovSave210(iUseRetrParam,iUseRetrParam);
  %%   r    = rSave210(iUseRetrParam,iUseRetrParam);
  %% elseif iSequential == 214
  %%   rcov = rcovSave214(iUseRetrParam,iUseRetrParam);
  %%   r    = rSave214(iUseRetrParam,iUseRetrParam);
  %% end
  %-------------------------
  
  % whos seSave rCovSave kSave raBTdeltanSave xbSave fuse iUseChan r k inv_se rcov rc
  chisqr0  = nansum(raBTdeltan'.*raBTdeltan');          %% only ~100 iSeqential chans
  xchisqr0 = nansum(raBTdeltanSave'.*raBTdeltanSave');  %% strow ~500 chans

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
      dx2 = k' * inv_se * raBTdeltan - r*(xn-xb);
    elseif invtype == 3
      dx2 = k'/inv_seF * raBTdeltan - r*(xn-xb);
    end
    deltax = dx1*dx2;
    figure(4); plot(diag(dx1)); colorbar;                              title('log10(dx1)');
    figure(5); imagesc(log10(abs(dx1))); colorbar;                     title('dx1');
    figure(6); plot(dx2);                                              title('dx2');
    figure(7); plot(fuse(iUseChan),raBTdeltan); plotaxis2;             title('raBTdeltaBT to fit')
    figure(8); plot(deltax.*driver.qrenorm(iUseRetrParam)'); plotaxis2; grid minor;   title('deltax.*qrenorm')
  
    if iDebug > 0

      %addpath /home/sergio/MATLABCODE; 
      figure(1); plot(f(inds(iUseChan)),k); grid
      figure(1); plot(f(inds(iUseChan)),k(:,1:5)); grid
      figure(2); pcolor(inv_se); shading flat; colorbar
      figure(2); plot(1:length(raBTdeltan),1./sqrt(diag(inv_se)),'b',1:length(raBTdeltan),-1./sqrt(diag(inv_se)),'b',1:length(raBTdeltan),raBTdeltan,'r'); grid
      figure(2); plot(1:length(raBTdeltan),1./sqrt(diag(inv_se)),'b',1:length(raBTdeltan),-1./sqrt(diag(inv_se)),'b'); grid
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

    if iSequential == 150
      if iFitTraceGas == -1
        fprintf(1,'deltax(SST) = %8.5f \n', deltax(1).*qrenorm(1))
      else
        fprintf(1,'deltax(CO2) = %8.5f \n', deltax(1).*qrenorm(1))
        fprintf(1,'deltax(SST) = %8.5f \n', deltax(2).*qrenorm(2))
      end
    end

    rodgers_rate = real(xn + deltax);
    figure(9); plot(1:length(xn),xn,'ko-',1:length(xn),deltax,'bx-',1:length(xn),real(xn+deltax),'r.-')
    figure(9); plot(1:length(xn),xn.*driver.qrenorm(iUseRetrParam)','ko-',1:length(xn),deltax.*driver.qrenorm(iUseRetrParam)','bx-',1:length(xn),real(xn+deltax).*driver.qrenorm(iUseRetrParam)','r.-')
      plotaxis2;   xlim([0 max(driver.jacobian.scalar_i)+1])
    hl = legend('orig xn','delta xn','new xn = (orig+delta)','location','best','fontsize',8); 
  
    ah0 = xn.*driver.qrenorm(iUseRetrParam)';
    dah = deltax.*driver.qrenorm(iUseRetrParam)';
    ah1 = real(xn+deltax).*driver.qrenorm(iUseRetrParam)';

    figure(12);
    if ~isfield(aux,'pavg')
      if iSequential == -1
        subplot(131); plot(ah0(driver.jacobian.water_i),1:length(driver.jacobian.water_i),'b',ah1(driver.jacobian.water_i),1:length(driver.jacobian.water_i),'r.-'); title('WV'); plotaxis2; set(gca,'ydir','reverse')
          ylim([1 length(driver.jacobian.water_i)+1]);
        subplot(132); plot(ah0(driver.jacobian.temp_i),1:length(driver.jacobian.water_i),'b',ah1(driver.jacobian.temp_i),1:length(driver.jacobian.water_i),'r.-');   title('T');  plotaxis2; set(gca,'ydir','reverse')
          ylim([1 length(driver.jacobian.water_i)+1]);
        subplot(133); plot(ah0(driver.jacobian.ozone_i),1:length(driver.jacobian.water_i),'b',ah1(driver.jacobian.ozone_i),1:length(driver.jacobian.water_i),'r.-'); title('O3'); plotaxis2; set(gca,'ydir','reverse')
          ylim([1 length(driver.jacobian.water_i)+1]);
      elseif iSequential == 60
        subplot(131); plot(ah0,1:length(driver.jacobian.water_i),'b',ah1,1:length(driver.jacobian.water_i),'r.-'); title('WV'); plotaxis2; set(gca,'ydir','reverse')
          ylim([1 length(driver.jacobian.water_i)+1]);
      elseif iSequential == 150
        if iFitTraceGas == -1
          subplot(132); plot(ah0(1),length(ah0),'bo',ah0(2:length(ah0)),1:length(ah0)-1,'b',ah1(1),length(ah0),'ro',ah1(2:length(ah0)),1:length(ah0)-1,'r.-');  title('T'); plotaxis2; set(gca,'ydir','reverse')
            ylim([1 length(driver.jacobian.water_i)+1]);
        else
          subplot(132); plot(ah0(2),length(ah0)-1,'bo',ah0(3:length(ah0)),1:length(ah0)-2,'b',ah1(2),length(ah0)-1,'ro',ah1(3:length(ah0)),1:length(ah0)-2,'r.-');  title('T'); plotaxis2; set(gca,'ydir','reverse')
            ylim([1 length(driver.jacobian.water_i)+1]);
        end
      elseif iSequential == 210
        junk = 2:3;
        subplot(132); plot(ah0(1),length(driver.jacobian.water_i)+1,'bo',ah0(6:length(ah0)),1:length(ah0)-5,'b',...
                           ah1(1),length(driver.jacobian.water_i)+1,'ro',ah1(6:length(ah0)),1:length(ah0)-5,'r.-');  title('T'); plotaxis2; set(gca,'ydir','reverse')
          ylim([1 length(driver.jacobian.water_i)+1]);
      elseif iSequential == 214
        junk = 2:iNXYZLay+1;
        subplot(131); plot(ah0(junk),1:iNXYZLay,'b',ah1(junk),1:iNXYZLay,'r.-');       title('WV'); plotaxis2; set(gca,'ydir','reverse'); ylim([1 iNXYZLay+1]);
        junk = (2:iNXYZLay+1)+iNXYZLay;
        subplot(132); plot(ah0(1),iNXYZLay+1,'bo',ah0(junk),1:iNXYZLay,'b',ah1(1),iNXYZLay+1,'ro',ah1(junk),1:iNXYZLay,'r.-');  title('T'); plotaxis2; set(gca,'ydir','reverse')
      elseif iSequential == 100
        subplot(133); plot(ah0,1:length(driver.jacobian.ozone_i),'b',ah1,1:length(driver.jacobian.ozone_i),'r.-'); title('O3'); plotaxis2; set(gca,'ydir','reverse')
          ylim([1 length(driver.jacobian.water_i)+1]);
      end

    elseif isfield(aux,'pavg')
      if iSequential == -1
        subplot(131); plot(ah0(driver.jacobian.water_i),aux.pavg,'b',ah1(driver.jacobian.water_i),aux.pavg,'r.-'); title('WV'); plotaxis2; set(gca,'ydir','reverse')
          ylim([1 aux.spres]);
        subplot(132); plot(ah0(driver.jacobian.temp_i),aux.pavg,'b',ah1(driver.jacobian.temp_i),aux.pavg,'r.-');   title('T');  plotaxis2; set(gca,'ydir','reverse')
          ylim([1 aux.spres]);
        subplot(133); plot(ah0(driver.jacobian.ozone_i),aux.pavg,'b',ah1(driver.jacobian.ozone_i),aux.pavg,'r.-'); title('O3'); plotaxis2; set(gca,'ydir','reverse')
          ylim([1 aux.spres]);
      elseif iSequential == 60
        subplot(131); plot(ah0,aux.pavg,'b',ah1,aux.pavg,'r.-'); title('WV'); plotaxis2; set(gca,'ydir','reverse')
          ylim([1 aux.spres]);
      elseif iSequential == 150
        if iFitTraceGas == -1
          subplot(132); plot(ah0(1),aux.spres,'bo',ah0(2:length(ah0)),aux.pavg,'b',ah1(1),aux.spres,'ro',ah1(2:length(ah0)),aux.pavg,'r.-');  title('T'); plotaxis2; set(gca,'ydir','reverse')
            ylim([1 aux.spres]);
        else
          subplot(132); plot(ah0(2),aux.spres,'bo',ah0(3:length(ah0)),aux.pavg,'b',ah1(2),aux.spres,'ro',ah1(3:length(ah0)),aux.pavg,'r.-');  title('T'); plotaxis2; set(gca,'ydir','reverse')
            ylim([1 aux.spres]);
        end
      elseif iSequential == 210
        junk = 2:3;
        subplot(132); plot(ah0(1),length(driver.jacobian.water_i)+1,'bo',ah0(6:length(ah0)),1:length(ah0)-5,'b',...
                           ah1(1),length(driver.jacobian.water_i)+1,'ro',ah1(6:length(ah0)),1:length(ah0)-5,'r.-');  title('T'); plotaxis2; set(gca,'ydir','reverse')
          ylim([1 aux.spres]);
      elseif iSequential == 214
        junk = 2:iNXYZLay+1;
        subplot(131); plot(ah0(junk),1:iNXYZLay,'b',ah1(junk),1:iNXYZLay,'r.-');       title('WV'); plotaxis2; set(gca,'ydir','reverse'); ylim([1 iNXYZLay+1]);
        junk = (2:iNXYZLay+1)+iNXYZLay;
        subplot(132); plot(ah0(1),iNXYZLay+1,'bo',ah0(junk),1:iNXYZLay,'b',ah1(1),iNXYZLay+1,'ro',ah1(junk),1:iNXYZLay,'r.-');  title('T'); plotaxis2; set(gca,'ydir','reverse')
      elseif iSequential == 100
        subplot(133); plot(ah0,aux.pavg,'b',ah1,aux.pavg,'r.-'); title('O3'); plotaxis2; set(gca,'ydir','reverse')
          ylim([1 aux.spres]);
      end
    end
  
    figure(1); 
    xn = rodgers_rate;
    xsave(iiS,ii,iUseRetrParam) = rodgers_rate;  %%% <<<< save rodgers_rate at iteration "for iiS = 1 : length(iaSequential)" "for ii = 1 : driver.oem.nloop"  >>>
  
    if ii <= driver.oem.nloop
      %% so this will be executed even if driver.oem.nloop == 1
  
      raBTdeltanIN = raBTdeltan;
      xn = rodgers_rate;
  
      % Form the computed rates; also see lines 108-121 of oem_lls.m
      thefitr      = zeros(1,length(driver.rateset.rates)); %% this uses all params including CO2/N2O/CH4; this uses all params including CO2/N2O/CH4; all 2645 chans
      thefitrdelta = zeros(1,length(driver.rateset.rates)); %% this just uses the params used in this iSequential round
                                                            %% all 2645 chans DOES NOT USE EG trace gas sometimes, or T(z)/ST sometimes, or WV sometimes , or O3 sometimes

      yn = xbSave;
      yn(iUseRetrParam) = xn;            %%% this is needed for next iSequential
      xsave(iiS,ii,:) = yn;              %%% <<<< save rodgers_rate and chisqr at iteration ii,iSequential of driver.oem.loop >>>>

      for ix = 1 : length(xbSave)
        %% this uses all params including CO2/N2O/CH4; this uses all params including CO2/N2O/CH4; all 2645 chans
        thefitr = thefitr + yn(ix)*m_ts_jac(:,ix)';
      end
      for ix = 1 : length(deltax)
        %% this just uses the params used in this iSequential round, all 2645 chans DOES NOT USE EG trace gas sometimes, or T(z)/ST sometimes, or WV sometimes , or O3 sometimes
        thefitrdelta = thefitrdelta + deltax(ix)*m_ts_jac(:,iUseRetrParam(ix))';
      end
      figure(7); plot(fuse(iUseChan),raBTdeltan,'cx-',fuse(iUseChan),thefitrdelta(inds(iUseChan)),'k',fuse(iUseChan),raBTdeltan'-thefitrdelta(inds(iUseChan)),'r','linewidth',2); plotaxis2; 
      title(['iSequential = ' num2str(iSequential) ' \newline (c) raBTdeltaN to be fitted (k) fit (r) diff'])
      hold on
        plot(f(inds),+driver.rateset.unc_rates(inds),'color',[1 1 1]*0.75);
        plot(f(inds),-driver.rateset.unc_rates(inds),'color',[1 1 1]*0.75);
        ylim([-1 +1]*0.05)
      hold off
    
      % Compute chisqr, and new raBTdeltan
      raBTdeltanIn   = raBTdeltan;                                              %% iSequential ~100 channels

      raBTdeltan     = driver.rateset.rates - thefitr';                         %% till  Jan 2021, 2645 chans
      raBTdeltanSave = raBTdeltan(inds);                                        %% being silly, but works beter? nah cant be, strow ~500 chans

      raBTzeltan = (driver.rateset.rates - tracegas_offset00) - thefitrdelta';  %% after Feb 2021, 2645 chans, why take out on;y what you have fitted huh?? silly
      raBTzeltan = (driver.rateset.rates - 1*tracegas_offset00) - thefitr';     %% after Jan 2023, remember we are assuming CO2/N2O/CH4 is fixed and all problems are with T/WV so should we take this out??? THAT IS UGH silly
      raBTzeltan = (driver.rateset.rates - 0*tracegas_offset00) - thefitr';     %% after Jan 2023, remember we are assuming CO2/N2O/CH4 is fixed and all problems are with T/WV so should we take this out??? NO, AND IDENTICAL to raBTdeltan
      raBTdeltanSave = raBTzeltan(inds);                                        %% see eg common_rodgers_initializations1.m, note raBTzeltan is different here than in rodgers.m, 

      raBTdeltaIterate(:,iiS+1) = raBTdeltanSave;                               %% strow ~500 chans
      xchisqr(iiS,ii)           = nansum(raBTdeltanSave'.*raBTdeltanSave');     %% strow ~500 chans

      raBTdeltan     = raBTdeltan(inds);                                        %% iSequential ~100 chans
      raBTzeltan     = raBTzeltan(inds);                                        %% iSequential ~100 chans
      chisqr(iiS,ii) = nansum(raBTzeltan'.*raBTzeltan');                        %% iSequential ~100 chans
     
  figure(22);
    plot(fuse,driver.rateset.rates(inds),'b-',fuse(iUseChan),driver.rateset.rates(inds(iUseChan)),'k.-',fuse(iUseChan),raBTdeltanIn,'g',fuse,raBTdeltanSave,'r')
      plotaxis2; hl = legend('all Strow chans rates','Sequential chans rates','starting raBTdelta','ending raBTdelta','location','best','fontsize',8);
  figure(23);
    if iSequential == -1
      junkST  = k(:,6);
      junkWV  = k(:,driver.jacobian.water_i);
      junkTz  = k(:,driver.jacobian.temp_i);
      junkO3  = k(:,driver.jacobian.ozone_i);
      plot(fuse(iUseChan),junkST,fuse(iUseChan),sum(junkWV,2),fuse(iUseChan),sum(junkTz,2),fuse(iUseChan),sum(junkO3,2)); plotaxis2; hl = legend('ST','WV','T','O3','location','best','fontsize',10); title('JAC')

      figure(2); clf; plot(deltax(driver.jacobian.water_i).*driver.qrenorm(driver.jacobian.water_i)',1:length(driver.jacobian.water_i),'b',...
                           yn(driver.jacobian.water_i).*driver.qrenorm(driver.jacobian.water_i)',1:length(driver.jacobian.water_i),'r');   plotaxis2; set(gca,'ydir','reverse'); title('\delta WV and current WV');
      figure(3); clf; plot(deltax(driver.jacobian.temp_i).*driver.qrenorm(driver.jacobian.temp_i)',1:length(driver.jacobian.temp_i),'b',...
                           yn(driver.jacobian.temp_i).*driver.qrenorm(driver.jacobian.temp_i)',1:length(driver.jacobian.temp_i),'r');      plotaxis2; set(gca,'ydir','reverse'); title('\delta Tz and current T');
      figure(4); clf; plot(deltax(driver.jacobian.ozone_i).*driver.qrenorm(driver.jacobian.ozone_i)',1:length(driver.jacobian.ozone_i),'b',...
                           yn(driver.jacobian.ozone_i).*driver.qrenorm(driver.jacobian.ozone_i)',1:length(driver.jacobian.ozone_i),'r');  plotaxis2; set(gca,'ydir','reverse'); title('\delta O3 and current O3');

    elseif iSequential == 150
      if iFitTraceGas == -1
        junkST  = k(:,1);
        junkTz  = k(:,2:length(driver.jacobian.temp_i)+1);
        plot(fuse(iUseChan),junkST,fuse(iUseChan),sum(junkTz,2)); plotaxis2; hl = legend('ST','T','location','best','fontsize',10); title('JAC')
        figure(3); clf; plot(deltax(2:length(deltax)).*driver.qrenorm(driver.jacobian.temp_i)',1:length(driver.jacobian.temp_i),'b',...
                             yn(driver.jacobian.temp_i).*driver.qrenorm(driver.jacobian.temp_i)',1:length(driver.jacobian.temp_i),'r');     plotaxis2; set(gca,'ydir','reverse'); title('\delta Tz');
      else
        junkST  = k(:,2);
        junkTz  = k(:,3:length(driver.jacobian.temp_i)+2);
        plot(fuse(iUseChan),junkST,fuse(iUseChan),sum(junkTz,2)); plotaxis2; hl = legend('ST','T','location','best','fontsize',10); title('JAC')
        figure(3); clf; plot(deltax(3:length(deltax)).*driver.qrenorm(driver.jacobian.temp_i)',1:length(driver.jacobian.temp_i),'b',...
                             yn(driver.jacobian.temp_i).*driver.qrenorm(driver.jacobian.temp_i)',1:length(driver.jacobian.temp_i),'r');     plotaxis2; set(gca,'ydir','reverse'); title('\delta Tz');
      end
    elseif iSequential == 100
      junkO3  = k;
      plot(fuse(iUseChan),sum(junkO3,2)); plotaxis2; hl = legend('O3','location','best','fontsize',10); title('JAC')
      figure(4); clf; plot(deltax.*driver.qrenorm(driver.jacobian.ozone_i)',1:length(driver.jacobian.ozone_i),'b',...
                           yn(driver.jacobian.ozone_i).*driver.qrenorm(driver.jacobian.ozone_i)',1:length(driver.jacobian.ozone_i),'r'); plotaxis2; set(gca,'ydir','reverse'); title('\delta O3');
    elseif iSequential == 60
      junkWV  = k;
      plot(fuse(iUseChan),sum(junkWV,2)); plotaxis2; hl = legend('WV','location','best','fontsize',10); title('JAC')
      figure(2); clf; plot(deltax.*driver.qrenorm(driver.jacobian.water_i)',1:length(driver.jacobian.water_i),'b',...
                           yn(driver.jacobian.water_i).*driver.qrenorm(driver.jacobian.water_i)',1:length(driver.jacobian.water_i),'r');   plotaxis2; set(gca,'ydir','reverse'); title('\delta WV');
    end
    figure(24); plot(fuse,raBTdeltaIterate,'linewidth',2); plotaxis2; hl = legend(num2str([0 [iaSequential(1:iiS)]]'),'location','best','fontsize',8); title('delta(obs-cal)')
    pause(0.1)

    iYesPlot = 1;    
    if driver.oem.doplots > 0 | iYesPlot > 0
      figure(10); clf
      figure(10); plot(f(inds),driver.rateset.rates(inds),'b.-',f(inds),driver.rateset.rates(inds) - tracegas_offset00(inds),'c.-',f(inds),thefitrdelta(inds),'k.-',f(inds),raBTzeltan,'r.-','linewidth',2); plotaxis2;
        hl = legend('input rates','signal=to fit after subtracting trace gas jacs','fit','signal''-fit','location','best','fontsize',8);
        title(['rodgers.m : ADJ SIGNAL \newline obs - fit at iteration ' num2str(ii)]); pause(0.1)
  
      figure(11); clf
      figure(11); plot(f(inds),driver.rateset.rates(inds),'b.-',f(inds),thefitr(inds),'k.-',f(inds),raBTdeltan,'r.-','linewidth',2); plotaxis2;
        hl = legend('input rates','fit','signal-fit','location','best','fontsize',8);
        title(['rodgers.m : RAW SIGNAL \newline obs - fit at iteration ' num2str(ii)]); pause(0.1)
      %plot(f(inds),raBTdeltanIN,f(inds),raBTdeltan,'r'); 
      %title(['obs - fit at iteration ' num2str(ii)]); pause(0.1)
    end
    xnIN = xn;
    end
  end

  fprintf(1,'iSequential = %3i : GND-2/GN-1/GND   WV and T, ST before \n',iSequential);
    junk = [xbSave(driver.jacobian.water_i(end-3):driver.jacobian.water_i(end)); xbSave(driver.jacobian.temp_i(end-3):driver.jacobian.temp_i(end)); xbSave(6)];
    fprintf(1,'%8.6f %8.6f %8.6f %8.6f      %8.6f %8.6f %8.6f %8.6f     %8.6f \n',junk);
  xbSave(iUseRetrParam) = xn; %% this is needed for next iSequential
  xnSave(iUseRetrParam) = xn; %% this is needed for next iSequential
  fprintf(1,'iSequential = %3i : GND-2/GN-1/GND   WV and T, ST after \n',iSequential);
    junk = [xbSave(driver.jacobian.water_i(end-3):driver.jacobian.water_i(end)); xbSave(driver.jacobian.temp_i(end-3):driver.jacobian.temp_i(end)); xbSave(6)];
    fprintf(1,'%8.6f %8.6f %8.6f %8.6f      %8.6f %8.6f %8.6f %8.6f     %8.6f \n',junk);

  figure(10); clf
  figure(10); plot(f(inds),driver.rateset.rates(inds),'b.-',f(inds),driver.rateset.rates(inds) - tracegas_offset00(inds),'c.-',f(inds),thefitrdelta(inds),'k.-',f(inds),raBTzeltan,'r.-','linewidth',2); plotaxis2;
          hl = legend('input rates','signal''= to fit after subtracting trace gas jacs','fit','signal''-fit','location','best','fontsize',8);
        title(['rodgers.m : ADJ SIGNAL \newline obs - fit at iteration ' num2str(ii)]); pause(0.1)
  
  figure(11); clf
  figure(11); plot(f(inds),driver.rateset.rates(inds),'b.-',f(inds),thefitr(inds),'k.-',f(inds),raBTdeltan,'r.-','linewidth',2); plotaxis2;
          hl = legend('input rates','fit','signal-fit','location','best','fontsize',8);
          title(['rodgers.m : RAW SIGNAL \newline obs - fit at iteration ' num2str(ii)]); pause(0.1)

  figure(11); title(['iSequential = ' num2str(iSequential)]); disp('rodgers_sequential.m 2'); ylim([-1 +1]*0.025); pause(0.1); 
    disp('exiting iSequential loop physical vars : [xbSave(ImportantParam) xnSave(ImportantParam) xb(ImportantParamX) xn(ImportantParamX)]')
    junk = [xbSave(iImportantParam).*driver.qrenorm(iImportantParam)' xnSave(iImportantParam).*driver.qrenorm(iImportantParam)' ...
     xb(iImportantParamX).*driver.qrenorm(iUseRetrParam(iImportantParamX))' xn(iImportantParamX).*driver.qrenorm(iUseRetrParam(iImportantParamX))'];
    for iijunk = 1 : 10
      fprintf(1,'%8.5f   %8.5f      %8.5f   %8.5f \n',junk(iijunk,:));
      if iijunk == 4 | iijunk == 6 | iijunk == 8
        disp('-------------------------------------------------------------')
      end
    end

  disp(' .................... >>>>>>>>>>> -------------------- <<<<<<<<<<<<< .....................')

  if iSequential >= -1 & iYesThisIsFine > 0
    wadachisqr = xchisqr(iiS,:);
    bestloop = find(wadachisqr ==  min(wadachisqr),1);
    wadaxsave  = squeeze(xsave(iiS,bestloop,:));
    fprintf(1,'for iiS = %2i iSequential = %2i bestloop (lowest chisqr) occured at iteration %3i \n',iiS,iSequential,bestloop)
    rodgers_rate_iiS(iiS,:) = wadaxsave;
    rodgers_chisqr_iiS(iiS) = wadachisqr(bestloop);
  end

  if driver.oem.nloop > 1
    fprintf(1,'printing out successive chisqr values (upto N-1 th iterate) ... %8.5f %8.5f %8.5f \n',[chisqr0 xchisqr(iiS,ii-1) xchisqr(iiS,ii)])
  elseif driver.oem.nloop == 1
    fprintf(1,'printing out successive chisqr values (upto N-1 th iterate) ... %8.5f %8.5f %8.5f \n',[chisqr0 chisqr0            xchisqr(iiS,ii)])
  end

  if iKeyBoard > 0
    keyboard_nowindow    
  else
    pause(0.1);
  end

  if iSequential == -1 & iiS == length(iaSequential)   
    disp('computing DofF')
    do_the_dof_avg_kernel
  end
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(1,'iSequential = %3i : GND-2/GN-1/GND   WV and T, ST exiting OEM \n',+8888);
  junk = [xbSave(driver.jacobian.water_i(end-3):driver.jacobian.water_i(end)); xbSave(driver.jacobian.temp_i(end-3):driver.jacobian.temp_i(end)); xbSave(6)];
  fprintf(1,'%8.6f %8.6f %8.6f %8.6f      %8.6f %8.6f %8.6f %8.6f     %8.6f \n',junk);

%% now find final best
fprintf(1,'iaSeqential chisqr = %5i %10.6f \n',[iaSequential; rodgers_chisqr_iiS]);
bestloop = find(rodgers_chisqr_iiS == min(rodgers_chisqr_iiS));
rodgers_rate = rodgers_rate_iiS(bestloop,:)';

fprintf(1,'iSequential = %3i : GND-2/GN-1/GND   WV and T, ST exiting OEM (best chisqr at stage %3i) \n',+9999,bestloop);
  junk = [rodgers_rate(driver.jacobian.water_i(end-3):driver.jacobian.water_i(end)); rodgers_rate(driver.jacobian.temp_i(end-3):driver.jacobian.temp_i(end)); rodgers_rate(6)];
  fprintf(1,'%8.6f %8.6f %8.6f %8.6f      %8.6f %8.6f %8.6f %8.6f     %8.6f \n',junk);

figure(24); plot(fuse,raBTdeltaIterate,'linewidth',2); plotaxis2; hl = legend(num2str([0 [iaSequential(1:iiS)]]'),'location','best','fontsize',8); title('delta(obs-cal)')
figure(2); hold on; plot(rodgers_rate(driver.jacobian.water_i).*driver.qrenorm(driver.jacobian.water_i)',1:length(driver.jacobian.water_i),'gx-'); plotaxis2; set(gca,'ydir','reverse'); hold off; title('\delta WV final');
figure(3); hold on; plot(rodgers_rate(driver.jacobian.temp_i).*driver.qrenorm(driver.jacobian.temp_i)',1:length(driver.jacobian.temp_i),'gx-');    plotaxis2; set(gca,'ydir','reverse'); hold off; title('\delta Tz final');
figure(4); hold on; plot(rodgers_rate(driver.jacobian.ozone_i).*driver.qrenorm(driver.jacobian.ozone_i)',1:length(driver.jacobian.ozone_i),'gx-'); plotaxis2; set(gca,'ydir','reverse'); hold off; title('\delta O3 final');

[mmmmS,nnnnS] = size(raBTdeltaIterate);
figure(25); plot(fuse,raBTdeltaIterate(:,nnnnS-1:nnnnS),'linewidth',2); plotaxis2; hl = legend(num2str([iaSequential(iiS-1:iiS)]'),'location','best','fontsize',8); title('delta(obs-cal)')

pause(0.1);
if iKeyBoard > 0
  ahF = real(rodgers_rate).*driver.qrenorm(iUseRetrParam)';
  figure(12); 
    if ~isfield(aux,'pavg')
      subplot(131); hold on; plot(ahF(driver.jacobian.water_i),1:length(driver.jacobian.water_i),'k.-'); title('WV'); plotaxis2; set(gca,'ydir','reverse'); hold off
      subplot(132); hold on; plot(ahF(driver.jacobian.temp_i),1:length(driver.jacobian.water_i),'k.-');   title('T'); plotaxis2; set(gca,'ydir','reverse'); hold off
      subplot(133); hold on; plot(ahF(driver.jacobian.ozone_i),1:length(driver.jacobian.water_i),'k.-'); title('O3'); plotaxis2; set(gca,'ydir','reverse'); hold off
    else
      subplot(131); hold on; plot(ahF(driver.jacobian.water_i),aux.pavg,'k.-'); title('WV'); plotaxis2; set(gca,'ydir','reverse'); hold off
      subplot(132); hold on; plot(ahF(driver.jacobian.temp_i),aux.pavg,'k.-');   title('T'); plotaxis2; set(gca,'ydir','reverse'); hold off
      subplot(133); hold on; plot(ahF(driver.jacobian.ozone_i),aux.pavg,'k.-'); title('O3'); plotaxis2; set(gca,'ydir','reverse'); hold off
    end
  keyboardstr = 'this is right before we exit rodgers_sequential.m, so check carefully'; keyboard_nowindow
end

%disp('lkjhsflKSHglshgslkhglhgslhgs')
%keyboard_nowindow
wxrodgers_rate = rodgers_rate;
rodgers_rate = rodgers_rate';   %% need this lousy line!!!!
fprintf(1,'iSequential = %3i : GND-2/GN-1/GND   WV and T, ST exiting OEM exiting OEM exiting OEM rodgers_sequential.m \n',+9999);
  junk = [wxrodgers_rate(driver.jacobian.water_i(end-3):driver.jacobian.water_i(end)); wxrodgers_rate(driver.jacobian.temp_i(end-3):driver.jacobian.temp_i(end)); wxrodgers_rate(6)];
  fprintf(1,'%8.6f %8.6f %8.6f %8.6f      %8.6f %8.6f %8.6f %8.6f     %8.6f \n',junk);

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
  
