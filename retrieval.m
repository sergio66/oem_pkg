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

  % are all the jacobians (B(T+dT)-B(T))/dt or (B(Q(1+dQ))-B(Q))/log(1+dQ))
  do_exp = [];
  if isfield(aux,'linearVSexpJAC')
    do_exp = find(aux.linearVSexpJAC == -1);
    fprintf(1,'found %3i variables that were done with logarithms for jacs \n',length(do_exp));
    driver.oem.finalrates(do_exp) = exp(driver.oem.finalrates(do_exp)) - 1;  %% SCALE FACTOR by which you multiply Qo, Aug 2018 I included the -1
    driver.oem.finalsigs(do_exp)  = exp(driver.oem.finalsigs(do_exp)) - 1;   %% SCALE FACTOR by which you multiply Qo, Aug 2018 I included the -1
  end

  fprintf(1,'IN MATLABCODE/oem_pkg/retrieval.m GND-2/GN-1/GND   WV and T, ST      and the uncertainty, after calling oem_lls');
  junk = [driver.oem.finalrates(driver.jacobian.water_i(end-3):driver.jacobian.water_i(end)); driver.oem.finalrates(driver.jacobian.temp_i(end-3):driver.jacobian.temp_i(end)); driver.oem.finalrates(6)];
  fprintf(1,'%8.6f %8.6f %8.6f %8.6f      %8.6f %8.6f %8.6f %8.6f     %8.6f \n',junk);
  junk = [driver.oem.finalsigs(driver.jacobian.water_i(end-3):driver.jacobian.water_i(end)); driver.oem.finalsigs(driver.jacobian.temp_i(end-3):driver.jacobian.temp_i(end)); driver.oem.finalsigs(6)];
  fprintf(1,'%8.6f %8.6f %8.6f %8.6f      %8.6f %8.6f %8.6f %8.6f     %8.6f \n',junk);
   
  % Get rid of OEM un-normalized coefficients
  driver.oem = rmfield(driver.oem,'coeffs');
  driver.oem = rmfield(driver.oem,'coeffssig');
end
