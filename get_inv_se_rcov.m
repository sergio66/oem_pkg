%{
 inv_se = inv(se);      x0 = norm(eye(size(se)) - inv_se * se,'fro');
 inv_se = pinv(se);     x1 = norm(eye(size(se)) - inv_se * se,'fro');
 %inv_se = invillco(se); x2 = norm(eye(size(se)) - inv_se * se,'fro');
 inv_se = inverse(se);  x3 = norm(eye(size(se)) - inv_se * se,'fro');
 %fprintf(1,'norms = %8.6f %8.6f %8.6f %8.6f \n',[x0 x1 x2 x3]);
 fprintf(1,'norms = %8.6f %8.6f %8.6f \n',[x0 x1    x3]);
 error('kjsf')
%}

% Do this once to save time, assume diagonal, no need for pinv   ORIG 201
if invtype == 0
  disp(' <<<<<<<< inv_se = inv(se)');            %%% NEW  pre Dec 2012
  inv_se = inv(se);         
elseif invtype == 1
  disp(' <<<<<<<< inv_se = pinv(se) DEFAULT');   %%% NEW  post Dec 2012
  inv_se = pinv(se);        
elseif invtype == 2
  disp(' <<<<<<<< inv_se = invillco(se)');       %%% NEW  post Apr 2019
  inv_se = invillco(se);    
elseif invtype == 3
  disp(' <<<<<<<< inv_se = inverse(se)');        %%% NEW  post Apr 2019
  inv_seF = factorize(se);  
  inv_se = inverse(se);     
  inv_se = inv_se * eye(size(inv_se));
elseif invtype == 4
  disp(' <<<<<<<< inv_se = inverse_ridge_regression_matrix(se)');            %%% NEW  pre Dec 2012
  inv_se = inverse_ridge_regression_matrix(se,kmax);   
elseif invtype == 5
  disp(' <<<<<<<< inv_se = inverse_minimum_eigenvalue_matrix(se)');            %%% NEW  pre Dec 2012
  kmaxrange   = 2 : 1 : 12;   kmaxrange = 10.^kmaxrange;
  sigminrange = -16 : 1 : -8; sigminrange = 10.^sigminrange;

  kmaxrange   = 2 : 0.25 : 12;   kmaxrange = 10.^kmaxrange;
  sigminrange = -16 : 0.25 : -8; sigminrange = 10.^sigminrange;

  kmaxrange   = 2 : 0.25 : 12;   kmaxrange = 10.^kmaxrange;
  sigminrange = -20 : 1 : -12; sigminrange = 10.^sigminrange;

  inv_se = inverse_minimum_eigenvalue_matrix_optim(se,kmaxrange,sigminrange,'Se');   
end
%if invtype <= 2
  oo = find(isinf(inv_se) | isnan(inv_se)); inv_se(oo) = 0;
%end

% driver.oem.cov needs to be inverted, since it is literally covariance; 
% the tikonov matrices below ("s") need not be inverted
if invtype == 0
  rcov = inv(driver.oem.cov);
elseif invtype == 1
  rcov = pinv(driver.oem.cov);
elseif invtype == 2
  rcov = invillco(driver.oem.cov);
elseif invtype == 3
  rcovF = factorize(driver.oem.cov);
  rcov  = inverse(driver.oem.cov);
  rcov  = rcov * eye(size(rcov));
elseif invtype == 4
  rcov = inverse_ridge_regression_matrix(driver.oem.cov,kmax);
elseif invtype == 5
  rcov = inverse_minimum_eigenvalue_matrix_optim(driver.oem.cov,kmaxrange,sigminrange,'driver.oem.cov');
end

% Use following line for only Tikhonov reg THIS SHOULD NOT BE INVERTED
l = get_l(driver.jacobian.numlays,1);    
s = transpose(l)*l;

%% now build the Tikhonov regularization block matrix, using "s"
lenS = length(driver.jacobian.scalar_i);
%% default : only column/stemp jacs, layer WV
rc = blkdiag(zeros(lenS,lenS),driver.oem.alpha_water*s);

%% always assumes you want to fit for WV ... may want to keep T fixed 9so no fit) and also may not want to fit O3
if isfield(driver.oem,'alpha_temp') & isfield(driver.oem,'alpha_ozone')
  rc = blkdiag(zeros(lenS,lenS),driver.oem.alpha_water*s,driver.oem.alpha_temp*s,driver.oem.alpha_ozone*s);
elseif isfield(driver.oem,'alpha_temp') & ~isfield(driver.oem,'alpha_ozone')
  rc = blkdiag(zeros(lenS,lenS),driver.oem.alpha_water*s,driver.oem.alpha_temp*s);
elseif ~isfield(driver.oem,'alpha_temp') & isfield(driver.oem,'alpha_ozone')
  rc = blkdiag(zeros(lenS,lenS),driver.oem.alpha_water*s,driver.oem.alpha_ozone*s);
end

% if invtype == 3, rcov could be a class rather than double
switch driver.oem.reg_type
  case 'reg_and_cov'
    r = rcov + rc;
  case 'reg'
    r = rc;
  case 'cov'
    r = rcov;
  otherwise
    disp('Incorrect choice driver.oem.reg_type')
end
