function driver = retrieval(driver,aux)

% Do the retrievals
driver = oem_lls(driver,aux);

% Translation to user units
renormalize = driver.qrenorm';

% Renormalize LLS coefficients
driver.lls.finalrates = driver.lls.coeffs.*renormalize;
sigs  = (driver.lls.coeffssig(:,2)-driver.lls.coeffssig(:,1))/2;
driver.lls.finalsigs = sigs.*renormalize;

if driver.oem.dofit
   % Renormalize OEM coefficients
   driver.oem.finalrates = driver.oem.coeffs'.*renormalize;
   driver.oem.finalsigs  = driver.oem.coeffssig'.*renormalize;

   % Get rid of OEM un-normalized coefficients
   driver.oem = rmfield(driver.oem,'coeffs');
   driver.oem = rmfield(driver.oem,'coeffssig');
end

% Get rid of LLS un-normalized coefficients
driver.lls = rmfield(driver.lls,'coeffs');
driver.lls = rmfield(driver.lls,'coeffssig');

   fprintf(1,'co2 error = %8.6f  \n',2.2*driver.oem.finalsigs(1));
   fprintf(1,'o3 error = %8.6f  \n',driver.oem.finalsigs(2)) ;
   fprintf(1,'N2O error = %8.6f  \n',driver.oem.finalsigs(3)) ;
   fprintf(1,'ch4 error = %8.6f  \n',driver.oem.finalsigs(4)); 
   fprintf(1,'cfc11 error = %8.6f  \n',driver.oem.finalsigs(5)); 
   fprintf(1,'sst error = %8.6f  \n', driver.oem.finalsigs(6)); 

