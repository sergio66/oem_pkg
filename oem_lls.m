function driver = oem_lls(driver,m_ts_jac);

% Make sure obs are finite!
inds = driver.jacobian.chanset;
obs  = real(driver.rateset.rates(inds)); 
inds = inds(isfinite(obs));
if length(inds) < length(driver.jacobian.chanset)
  disp('Some input rates have NaNs, now removed')
  driver.jacobian.chanset = inds;
end

% Make sure obs lie between MIN and MAX
obs  = real(driver.rateset.rates(inds)); 
junk = find(obs < driver.rateset.max & obs > driver.rateset.min);
inds = inds(junk);
driver.jacobian.chanset = inds;
clear obs junk

% Make sure channel noise levels are less than 2 K, new July 2013
% (this is sometimes not satisfied at the 650-700 cm-1 AIRS channels)
obs  = real(driver.rateset.rates(inds)); 
junk = find(abs(driver.rateset.unc_rates(inds)) <= 1.0);
inds = inds(junk);
driver.jacobian.chanset = inds;
clear obs junk

% see if we can find a data spike
if driver.rateset.despike > 0
  junkFF = driver.f(inds);
  junkA  = real(driver.rateset.rates(inds));
  plot(0.5*(junkFF(1:end-1)+junkFF(2:end)),diff(junkA));
  [junkC,idxC] = despike2(real(driver.rateset.rates(inds)),driver.rateset.despike);
  if length(idxC) > 0
    disp(['  <<<<<<< found ' num2str(length(idxC)) ' spikes, getting rid of them ...'])
    indsX = setdiff(1:length(inds),idxC);
    indsX = inds(indsX);
  else
    indsX = inds;
  end
  inds = indsX;
  driver.jacobian.chanset = inds;
  plot(junkFF,junkA,'bo-',driver.f(inds),real(driver.rateset.rates(inds)),'rx-');
end

[mm,nn] = size(m_ts_jac(inds,:));
if (length(inds) < nn & driver.oem.dofit)
  fprintf(1,'  length(inds) = %4i \n',length(inds));
  fprintf(1,'  size jac mm,nn = %4i %4i\n',mm,nn);
  disp('Error: More jacobians than channels!')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LLS

[coeffs,coeffssig] = regress(driver.rateset.rates(inds),m_ts_jac(inds,:));

thefit = zeros(length(driver.rateset.rates),1);
for ix = 1 : length(coeffs)
  thefit = thefit + coeffs(ix)*m_ts_jac(:,ix);
end

if driver.jacobian.qstYesOrNo(1) == 1 & strfind(driver.jacobian.qstnames{1},'CO2')
  fprintf(1,'quick lls co2 estimate = %8.6f ppmv/yr \n',coeffs(1)*driver.qrenorm(1))
else
  disp('driver.jacobian.co2 or this is not co2 = false ... not printing result!')
end
junk = length(driver.jacobian.qstYesOrNo);
if driver.jacobian.qstYesOrNo(junk) == 1 & strfind(driver.jacobian.qstnames{junk},'stemp')
  fprintf(1,'quick stemp estimate = %8.6f K/yr \n',coeffs(junk)*driver.qrenorm(junk))
else
  disp('driver.jacobian.stemp = false or not stemp ... not printing result!')
end

driver.lls.coeffs    = coeffs;
driver.lls.coeffssig = coeffssig;
driver.lls.fit       = thefit;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OEM

% Don't do OEM fit if dofit is false.  This is for LLS studies
if driver.oem.dofit

  % Get the observed rate errors, which ultimately are a combination of
  %  (a) get_rates.m : driver.rateset.unc_rates = real(squeeze(b_err_obs(iibin,:,2))');
  %  (b) nc_rates.m  : ncerrors = real(nc' .* driver.rateset.unc_rates); where nc is 1, or comes from a file
  ncerrors = nc_rates(driver);

  % xset are the apriori rates, in use units eg ppm/yr, K/yr
  xb = load(driver.oem.apriori_filename,'apriori');
  xb = xb.apriori(driver.jacindex,:);
  [mm,nn] = size(xb);
  if nn > 1
    xb = xb(:,driver.iibin);
  end
  xb = xb./driver.qrenorm';

  % Form structure needed by rodgers.m
  aux_stuff.m_ts_jac = m_ts_jac;
  aux_stuff.xb       = xb;
  aux_stuff.ncerrors = ncerrors;

  % Do the OEM retrieval
  if driver.oem.rates == +1;
    %% do rates
    [rodgers_rate,errorx,dofs,gain,kern,inds,r,se,inv_se,se_errors] = rodgers(driver,aux_stuff);
  elseif driver.oem.rates == -1;
    %% do regular spectra
    [rodgers_rate,errorx,dofs,gain,kern,inds,r,se,inv_se,se_errors] = rodgers_spectra(driver,aux_stuff);
  end

  driver.jacobian.chanset_used = inds;

  %% show the terms used in the Se error matrix
  driver.oem.spectral_errors     = se_errors.ncerrors;
  driver.oem.forwardmodel_errors = se_errors.fmerrors;
  driver.oem.inv_se              = inv_se;
  driver.oem.se                  = se;
  
  %% show the "r" matrix used for regularization
  driver.oem.regularizationmatrix = r;

  % Build the output structure
  driver.oem.gain  = gain; 
  driver.oem.ak    = kern;
  driver.oem.dofs  = dofs;
  
  coeffsr          = rodgers_rate;
  coeffssigr       = sqrt(diag(errorx)');    % AVT sigs should be square root of the covariance diagonal

  % Form the computed rates
  thefitr = zeros(1,length(driver.rateset.rates));
  for ix = 1 : length(coeffs)
    thefitr = thefitr + coeffsr(ix)*m_ts_jac(:,ix)';
  end

  % Compute chisqr
  fit_minus_obs = thefitr - driver.rateset.rates'; 
  fit_minus_obs = fit_minus_obs(:,inds);
  chisqrr = sum(fit_minus_obs'.*fit_minus_obs');

  % Compute variance of fitted coefficients
  coeff_var = xb - coeffsr'; 
  coeff_var_qst = sum(coeff_var(driver.jacobian.iqst).^2);
  if driver.jacobian.numlays > 0 & driver.jacobian.numQlays > 0
    for ii = 1 : driver.jacobian.numQlays
      junk = ['coeff_var_Q' num2str(ii) ' = sum(coeff_var(driver.jacobian.iQ' num2str(ii) ').^2);']; 
      eval(junk)
    end
  end
  if driver.jacobian.numlays > 0
    coeff_var_temp = sum(coeff_var(driver.jacobian.itemp).^2);
  end

  % More output
  driver.oem.coeffs          = coeffsr;
  driver.oem.coeffssig       = coeffssigr;
  driver.oem.fit             = thefitr;
  driver.oem.chisqr          = chisqrr;
  driver.oem.coeff_var_qst   = coeff_var_qst;
  for ii = 1 : driver.jacobian.numQlays
    junk = ['driver.oem.coeff_var_Q' num2str(ii) '   = coeff_var_Q' num2str(ii) ';']; eval(junk)
  end
  if driver.jacobian.numlays > 0
    driver.oem.coeff_var_temp  = coeff_var_temp;
  end

  if driver.jacobian.qstYesOrNo(1) == 1 & strfind(driver.jacobian.qstnames{1},'CO2')
    fprintf(1,'quick oem co2 estimate = %8.6f ppmv/yr \n',coeffsr(1)*driver.qrenorm(1))
  else
    disp('driver.jacobian.co2 = false or this is not CO2... not printing result!')
  end

  if driver.jacobian.qstYesOrNo(2) == 1 & strfind(driver.jacobian.qstnames{2},'O3')
    fprintf(1,'quick oem o3 estimate = %8.6f ppmv/yr \n',coeffsr(2)*driver.qrenorm(2))
  else
    disp('driver.jacobian.o3 = false or this is not CO2... not printing result!')
  end

  if driver.jacobian.qstYesOrNo(3) == 1 & strfind(driver.jacobian.qstnames{3},'N2O')
    fprintf(1,'quick oem n2o estimate = %8.6f ppmv/yr \n',coeffsr(3)*driver.qrenorm(3))
  else
    disp('driver.jacobian.n2o = false or this is not N2O... not printing result!')
  end

  if driver.jacobian.qstYesOrNo(4) == 1 & strfind(driver.jacobian.qstnames{4},'CH4')
    fprintf(1,'quick oem ch4 estimate = %8.6f ppmv/yr \n',coeffsr(4)*driver.qrenorm(4))
  else
    disp('driver.jacobian.ch4 = false or this is not CH4... not printing result!')
  end

  if driver.jacobian.qstYesOrNo(5) == 1 & strfind(driver.jacobian.qstnames{5},'CFC11')
    fprintf(1,'quick oem cfc11 estimate = %8.6f ppmv/yr \n',coeffsr(5)*driver.qrenorm(5))
  else
    disp('driver.jacobian.cfc11 = false or this is not CFC11... not printing result!')
  end

  junk = length(driver.jacobian.qstYesOrNo);
  if driver.jacobian.qstYesOrNo(junk) == 1 & strcmp(driver.jacobian.qstnames{junk},'stemp')
    fprintf(1,'quick oem sst estimate = %8.6f K/yr \n',coeffsr(junk)*driver.qrenorm(junk))
  else
    disp('driver.jacobian.stemp = false or not stemp ... not printing result!')
  end

end % of do/don't do  OEM fit
