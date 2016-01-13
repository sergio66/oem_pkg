function driver = retrieval(driver,aux)

% Use sq2 = 2:200 when NOT fitting CO2.  This needs to be programmed into
% driver.oem.
%sq2 = 2:200;
%% default
sq2 = 1:length(driver.jacobian.scalar_i) + length(driver.jacobian.water_i) + length(driver.jacobian.temp_i);
if isfield(driver.jacobian,'ozone_i')
  sq2 = 1:length(driver.jacobian.scalar_i) + length(driver.jacobian.water_i) + length(driver.jacobian.temp_i) + length(driver.jacobian.ozone_i);
end

% Translation to user units
renormalize = driver.qrenorm(sq2)';
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
