function driver = retrieval(driver,aux)

% Translation to user units
renormalize = driver.qrenorm';
driver = oem_lls(driver,aux);

% Do the LLS retrievals
if driver.lls.dofit
   driver.lls.finalrates = driver.lls.coeffs.*renormalize;
   sigs  = (driver.lls.coeffssig(:,2)-driver.lls.coeffssig(:,1))/2;
   driver.lls.finalsigs = sigs.*renormalize;
% Get rid of LLS un-normalized coefficients
   driver.lls = rmfield(driver.lls,'coeffs');
   driver.lls = rmfield(driver.lls,'coeffssig');
end

if driver.oem.dofit
% Renormalize OEM coefficients
   driver.oem.finalrates = driver.oem.coeffs'.*renormalize;
   driver.oem.finalsigs  = driver.oem.coeffssig'.*renormalize;

   % Get rid of OEM un-normalized coefficients
   driver.oem = rmfield(driver.oem,'coeffs');
   driver.oem = rmfield(driver.oem,'coeffssig');
end
