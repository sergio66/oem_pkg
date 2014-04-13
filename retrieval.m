function driver = retrieval(driver,aux)

% Do the retrievals
driver = oem_lls(driver,aux);

% Translation to user units
renormalize = driver.qrenorm';

% Renormalize LLS coefficients
driver.lls.finalrates = driver.lls.coeffs.*renormalize;
sigs  = (driver.lls.coeffssig(:,2)-driver.lls.coeffssig(:,1))/2;
driver.lls.finalsigs = sigs.*renormalize;
% Get rid of LLS un-normalized coefficients
driver.lls = rmfield(driver.lls,'coeffs');
driver.lls = rmfield(driver.lls,'coeffssig');

if driver.oem.dofit
   % Renormalize OEM coefficients
   driver.oem.finalrates = driver.oem.coeffs'.*renormalize;
   driver.oem.finalsigs  = driver.oem.coeffssig'.*renormalize;

   % Get rid of OEM un-normalized coefficients
   driver.oem = rmfield(driver.oem,'coeffs');
   driver.oem = rmfield(driver.oem,'coeffssig');
end
