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
  driver.oem.spectral_errors     = aux.ncerrors;
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

% Sergio, do we need this?
%   % Compute variance of fitted coefficients
%   coeff_var = aux.xb - coeffsr'; 
%   coeff_var_qst = sum(coeff_var(driver.jacobian.iqst).^2);
%   if driver.jacobian.numlays > 0 & driver.jacobian.numQlays > 0
%     for ii = 1 : driver.jacobian.numQlays
%       junk = ['coeff_var_Q' num2str(ii) ' = sum(coeff_var(driver.jacobian.iQ' num2str(ii) ').^2);']; 
%       eval(junk)
%     end
%   end
%   if driver.jacobian.numlays > 0
%     coeff_var_temp = sum(coeff_var(driver.jacobian.itemp).^2);
%   end
  driver.oem.coeffs          = coeffsr;
  driver.oem.coeffssig       = coeffssigr;
  driver.oem.fit             = thefitr;
  driver.oem.chisqr          = chisqrr;

% Sergio: needed?  
%   driver.oem.coeff_var_qst   = coeff_var_qst;
%   driver.oem.coeff_var_temp  = coeff_var_temp;

end % do oem fit
