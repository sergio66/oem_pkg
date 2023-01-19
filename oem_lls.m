function driver = oem_lls(driver,aux);

%LLS-------------------------------------------------------------------------
if driver.lls.dofit

   chanset = driver.jacobian.chanset;
   [coeffs,coeffssig] = regress(driver.rateset.rates(chanset),aux.m_ts_jac(chanset,:));
   thefit = zeros(length(driver.rateset.rates),1);
   for ix = 1 : length(coeffs)
      thefit = thefit + coeffs(ix)*aux.m_ts_jac(:,ix);
   end
   driver.lls.coeffs    = coeffs;
   driver.lls.coeffssig = coeffssig;
   driver.lls.fit       = thefit;
end

%OEM--------------------------------------------------------------------------
% Don't do OEM fit if dofit is false.  This is for LLS studies
if driver.oem.dofit

  if ~isfield(driver,'iaSequential')
    disp('in oem_lls : you must have git checkout oldbranch : setting driver.iaSequential = -1')
    driver.iaSequential = -1;
  end

  % Do the OEM retrieval
  if driver.iaSequential(1) == -1 & length(driver.iaSequential) == 1
    %% default ie do all geophysical parameters in one massive gulp
    [rodgers_rate,errorx,dofs,cdofs,gain,kern,r,se,inv_se,se_errors,kern_water,kern_temp,kern_ozone,bestloop,raBTdeltan00] = rodgers(driver,aux);
  else
    %% do geophysical params sequentially 
    [rodgers_rate,errorx,dofs,cdofs,gain,kern,r,se,inv_se,se_errors,kern_water,kern_temp,kern_ozone,bestloop,raBTdeltan00] = rodgers_sequential(driver,aux);
  end

  % Save terms used in the Se error matrix
  driver.oem.forwardmodel_errors = se_errors.fmerrors;
  driver.oem.inv_se              = inv_se;
  driver.oem.se                  = se;
  
  % what did we actually fit (= input signal(v) -  sum_i(jac(v,i)*tracegas_rate(i)))
  driver.oem.spectral_deltan00 = raBTdeltan00;

  % Build the output structure
  driver.oem.gain  = gain; 
  driver.oem.ak    = kern;
  driver.oem.dofs  = dofs;
  driver.oem.cdofs = cdofs;
  driver.oem.ak_water = kern_water; 
  driver.oem.ak_temp  = kern_temp;
  driver.oem.ak_ozone = kern_ozone;
  driver.oem.error_cov = errorx;
  driver.oem.bestloop  = bestloop;

  coeffsr          = rodgers_rate;
  coeffssigr       = sqrt(diag(errorx)');    % AVT sigs should be square root of the covariance diagonal

  % Form the computed rates
  thefitr = zeros(1,length(driver.rateset.rates));
  for ix = 1 : length(coeffsr)
    thefitr = thefitr + coeffsr(ix)*aux.m_ts_jac(:,ix)';
    %if length(intersect(ix,[7 27 47])) == 1
    %  %% separate out rtace gas, WV, T, O3
    %  disp('-----------')
    %end
    %fprintf(1,'oem_lls 1231 cm-1 chID 1520  ix,truecoeff/qrenorm,truejac*qrenorm,trueproduct : %3i %10.6f %10.6e %10.6e \n',ix,coeffsr(ix),aux.m_ts_jac(1520,ix),coeffsr(ix)*aux.m_ts_jac(1520,ix))
  end
  %woo = aux.m_ts_jac; whos woo thefitr

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
