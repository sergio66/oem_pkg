% Error analysis and diagnostics  

%lala1 = diag(inv(AKstuff)*AKstuff);
%lala2 = diag(pinv(AKstuff)*AKstuff);
%lala3 = diag(invillco(AKstuff)*AKstuff);
%lala4 = diag(inv(U)*inv(L) * AKstuff);
%figure(5); plot(1:length(lala1),lala1,1:length(lala1),lala2,1:length(lala1),lala3); error('klfkf')


if invtype == 0
  errorx = inv(k' * inv_se * k + r);     %% decided pinv is too unstable, Aug 2018
elseif invtype == 1
  errorx = pinv(k' * inv_se * k + r);    %% decided pinv is too unstable, Nov 2013, but not much difference
elseif invtype == 2
  errorx = invillco(k' * inv_se * k + r);%% decided pinv is too unstable, Nov 2013, but not much difference
elseif invtype == 3
  errorx = double(k' * inv_se * k + r);  %% decided pinv is too unstable, May 2019, but not much difference
  errorx = inverse(errorx);
  errorx = errorx * eye(size(errorx));
elseif invtype == 4
  errorx = inverse_ridge_regression_matrix(k' * inv_se * k + r,kmax);     %% decided pinv is too unstable, Aug 2018
elseif invtype == 5
  errorx = inverse_minimum_eigenvalue_matrix_optim(k' * inv_se * k + r,kmaxrange,sigminrange,'errorx');     %% decided pinv is too unstable, Aug 2018
elseif invtype == 9999
  AKstuff = (k' * inv_se * k + r);
  [L,U] = lu(AKstuff);
  errorx = inv(U)*inv(L);                %% trying LU
end

dofsx  = errorx * r; 
dofsx  = eye(size(dofsx)) - dofsx; 
dofs   = trace(dofsx);
cdofs  = diag(dofsx);                 %% so we can do cumulative d.of.f

% Gain is relative weight of first guess and observations
r_water  = r(driver.jacobian.water_i,driver.jacobian.water_i); 
r_temp   = r(driver.jacobian.temp_i,driver.jacobian.temp_i); 
if invtype == 0
  inv_r       = inv(r);
  inv_r_water = inv(r_water); 
  inv_r_temp  = inv(r_temp); 
elseif invtype == 1
  inv_r       = pinv(r);
  inv_r_water = pinv(r_water); 
  inv_r_temp  = pinv(r_temp); 
elseif invtype == 2
  inv_r       = invillco(r);
  inv_r_water = invillco(r_water); 
  inv_r_temp  = invillco(r_temp); 
elseif invtype == 3
  inv_r       = inverse(r);
  inv_r_water = inverse(r_water); 
  inv_r_temp  = inverse(r_temp); 
elseif invtype == 4
  inv_r       = inverse_ridge_regression_matrix(r,kmax);
  inv_r_water = inverse_ridge_regression_matrix(r_water,kmax); 
  inv_r_temp  = inverse_ridge_regression_matrix(r_temp,kmax); 
elseif invtype == 5
  inv_r       = inverse_minimum_eigenvalue_matrix_optim(r,kmaxrange,sigminrange,'inv_r');
  inv_r_water = inverse_minimum_eigenvalue_matrix_optim(r_water,kmaxrange,sigminrange,'inv_r_water'); 
  inv_r_temp  = inverse_minimum_eigenvalue_matrix_optim(r_temp,kmaxrange,sigminrange,'inv_r_temp'); 
end

% inv operator seems OK for this matrix; if problems go back to pinv
k_water    = k(:,driver.jacobian.water_i); 
k_temp     = k(:,driver.jacobian.temp_i); 
if invtype == 0
  gain       = inv_r *k' * inv(k * inv_r * k' + se);
  gain_water = inv_r_water*k_water'*inv(k_water*inv_r_water*k_water'+se); 
  gain_temp  = inv_r_temp*k_temp'*inv(k_temp*inv_r_temp*k_temp'+se); 
elseif invtype == 1
  gain       = inv_r *k' * pinv(k * inv_r * k' + se);
  gain_water = inv_r_water*k_water'*pinv(k_water*inv_r_water*k_water'+se); 
  gain_temp  = inv_r_temp*k_temp'*pinv(k_temp*inv_r_temp*k_temp'+se); 
elseif invtype == 2
  gain       = inv_r *k' * invillco(k * inv_r * k' + se);
  gain_water = inv_r_water*k_water'*invillco(k_water*inv_r_water*k_water'+se); 
  gain_temp  = inv_r_temp*k_temp'*invillco(k_temp*inv_r_temp*k_temp'+se); 
elseif invtype == 3
  junk   = double(k * inv_r * k' + se);
    gain = inv_r *k' * inverse(junk);
  junk = double(k_water*inv_r_water*k_water'+se);
    gain_water = inv_r_water*k_water'*inverse(junk);
  junk = double(k_temp*inv_r_temp*k_temp'+se);
    gain_temp  = inv_r_temp*k_temp'*inverse(junk);
elseif invtype == 4
  gain       = inv_r *k' * inverse_ridge_regression_matrix(k * inv_r * k' + se,kmax);
  gain_water = inv_r_water*k_water'*inverse_ridge_regression_matrix(k_water*inv_r_water*k_water'+se,kmax); 
  gain_temp  = inv_r_temp*k_temp'*inverse_ridge_regression_matrix(k_temp*inv_r_temp*k_temp'+se,kmax); 
elseif invtype == 5
  gain       = inv_r *k' * inverse_minimum_eigenvalue_matrix_optim(k * inv_r * k' + se,kmaxrange,sigminrange,'gain');
  gain_water = inv_r_water*k_water'*inverse_minimum_eigenvalue_matrix_optim(k_water*inv_r_water*k_water'+se,kmaxrange,sigminrange,'gain_water'); 
  gain_temp  = inv_r_temp*k_temp'*inverse_minimum_eigenvalue_matrix_optim(k_temp*inv_r_temp*k_temp'+se,kmaxrange,sigminrange,'gain_temp'); 
end

% Compute averaging kernel
ak = gain * k;   
ak_water = gain_water*k_water; 
ak_temp  = gain_temp*k_temp; 

if isfield(driver.oem,'alpha_ozone')
  r_ozone  = r(driver.jacobian.ozone_i,driver.jacobian.ozone_i); 
  k_ozone    = k(:,driver.jacobian.ozone_i);
  if invtype == 0 
    inv_r_ozone = inv(r_ozone);
    gain_ozone = inv_r_ozone*k_ozone'*inv(k_ozone*inv_r_ozone*k_ozone'+se);
  elseif invtype == 1 
    inv_r_ozone = pinv(r_ozone);
    gain_ozone = inv_r_ozone*k_ozone'*pinv(k_ozone*inv_r_ozone*k_ozone'+se);
  elseif invtype == 2
    inv_r_ozone = invillco(r_ozone);
    gain_ozone = inv_r_ozone*k_ozone'*invillco(k_ozone*inv_r_ozone*k_ozone'+se);
  elseif invtype == 3
    inv_r_ozone = inverse(r_ozone);
    gain_ozone = inv_r_ozone*k_ozone'*inverse(k_ozone*inv_r_ozone*k_ozone'+se);
  elseif invtype == 4
    inv_r_ozone = inverse_ridge_regression_matrix(r_ozone,kmax);
    gain_ozone = inv_r_ozone*k_ozone'*inverse_ridge_regression_matrix(k_ozone*inv_r_ozone*k_ozone'+se,kmax);
  elseif invtype == 5
    inv_r_ozone = inverse_minimum_eigenvalue_matrix_optim(r_ozone,kmaxrange,sigminrange,'inv_r_ozone');
    gain_ozone = inv_r_ozone*k_ozone'*inverse_minimum_eigenvalue_matrix_optim(k_ozone*inv_r_ozone*k_ozone'+se,kmaxrange,sigminrange,'gain_ozone');
  end
  ak_ozone = gain_ozone*k_ozone;
else
  ak_ozone = zeros(size(ak_water));
end
