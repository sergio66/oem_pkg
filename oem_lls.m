function driver = oem_lls(driver,m_ts_jac);

% Make sure obs are finite!
inds = driver.jacobian.chanset;
obs  = real(driver.rateset.rates(inds)); 
inds = inds( isfinite(obs) );
if length(inds) < length(driver.jacobian.chanset)
  disp('Some input rates have NaNs, now removed')
end
obs  = real(driver.rateset.rates(inds)); 
junk = find(obs < driver.rateset.max & obs > driver.rateset.min);
inds = inds(junk);
driver.jacobian.chanset = inds;
clear obs junk

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

  % Get the observed rate errors
  ncerrors = nc_rates(driver);

  % xset are the apriori rates, in use units eg ppm/yr, K/yr
  xb = load(driver.oem.apriori_filename,'apriori');
  xb = xb.apriori(driver.jacindex)./driver.qrenorm';

  % Form structure needed by rodgers.m
  aux_stuff.m_ts_jac = m_ts_jac;
  aux_stuff.xb       = xb;
  aux_stuff.ncerrors = ncerrors;

  % Do the OEM retrieval
  [rodgers_rate,errorx,dofs,gain,kern,inds] = rodgers(driver,aux_stuff);
  driver.jacobian.chanset_used = inds;

  % Build the output structure
  driver.oem.gain  = gain; 
  driver.oem.ak    = kern;
  driver.oem.dofs  = dofs;
  coeffsr          = rodgers_rate;
  coeffssigr       = diag(errorx)';

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
  for ii = 1 : driver.jacobian.numQlays
    junk = ['coeff_var_Q' num2str(ii) ' = sum(coeff_var(driver.jacobian.iQ' num2str(ii) ').^2);']; eval(junk)
  end
  coeff_var_temp = sum(coeff_var(driver.jacobian.itemp).^2);

  % More output
  driver.oem.coeffs          = coeffsr;
  driver.oem.coeffssig       = coeffssigr;
  driver.oem.fit             = thefitr;
  driver.oem.chisqr          = chisqrr;
  driver.oem.coeff_var_qst   = coeff_var_qst;
  for ii = 1 : driver.jacobian.numQlays
    junk = ['driver.oem.coeff_var_Q' num2str(ii) '   = coeff_var_Q' num2str(ii) ';']; eval(junk)
  end
  driver.oem.coeff_var_temp  = coeff_var_temp;

  if driver.jacobian.qstYesOrNo(1) == 1 & strfind(driver.jacobian.qstnames{1},'CO2')
    fprintf(1,'quick oem co2 estimate = %8.6f ppmv/yr \n',coeffsr(1)*driver.qrenorm(1))
  else
    disp('driver.jacobian.co2 = false or this is not CO2... not printing result!')
  end

  junk = length(driver.jacobian.qstYesOrNo);
  if driver.jacobian.qstYesOrNo(junk) == 1 & strcmp(driver.jacobian.qstnames{junk},'stemp')
    fprintf(1,'quick oem sst estimate = %8.6f K/yr \n',coeffsr(junk)*driver.qrenorm(junk))
  else
    disp('driver.jacobian.stemp = false or not stemp ... not printing result!')
  end

end % of do/don't do  OEM fit