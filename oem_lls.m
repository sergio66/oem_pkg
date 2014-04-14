function driver = oem_lls(driver,aux);

%LLS-------------------------------------------------------------------------
[coeffs,coeffssig] = regress(driver.rateset.rates,aux.m_ts_jac);
thefit = zeros(length(driver.rateset.rates),1);
for ix = 1 : length(coeffs)
  thefit = thefit + coeffs(ix)*aux.m_ts_jac(:,ix);
end
driver.lls.coeffs    = coeffs;
driver.lls.coeffssig = coeffssig;
driver.lls.fit       = thefit;

%OEM--------------------------------------------------------------------------
% Don't do OEM fit if dofit is false.  This is for LLS studies
if driver.oem.dofit

  % Do the OEM retrieval
  [rodgers_rate,errorx,dofs,cdofs,gain,kern,r,se,inv_se,se_errors,kern_water,kern_temp] = rodgers(driver,aux);

  % Save terms used in the Se error matrix
  driver.oem.spectral_errors     = driver.rateset.unc_rates;
  driver.oem.forwardmodel_errors = se_errors.fmerrors;
  driver.oem.inv_se              = inv_se;
  driver.oem.se                  = se;
  
  % Save actual r matrix
  driver.oem.regularizationmatrix = r;

  % Build the output structure
  driver.oem.gain  = gain; 
  driver.oem.ak    = kern;
  driver.oem.dofs  = dofs;
  driver.oem.cdofs = cdofs;
  driver.oem.ak_water = kern_water; 
  driver.oem.ak_temp = kern_temp; 
  coeffsr          = rodgers_rate;
  coeffssigr       = sqrt(diag(errorx)');    % AVT sigs should be square root of the covariance diagonal

  % Form the computed rates
  thefitr = zeros(1,length(driver.rateset.rates));
  for ix = 1 : length(coeffs)
    thefitr = thefitr + coeffsr(ix)*aux.m_ts_jac(:,ix)';
  end

  inds     = driver.jacobian.chanset;
  % Compute chisqr
  fit_minus_obs = thefitr - driver.rateset.rates'; 
  fit_minus_obs = fit_minus_obs(:,inds);
  chisqrr = sum(fit_minus_obs'.*fit_minus_obs');
  driver.oem.coeffs          = coeffsr;
  driver.oem.coeffssig       = coeffssigr;
  driver.oem.fit             = thefitr;
  driver.oem.chisqr          = chisqrr;

end % do oem fit
