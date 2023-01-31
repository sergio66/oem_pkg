if iSequential == -1
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif iSequential == -1
  X2 = driver.oem.cov;
  reducefact = 1;
  reducefact = 10;
  reducefact = sqrt(10);
  X2 = X2/(reducefact*reducefact); %% make unc^2 in parameters reducefact times 

  % driver.oem.cov needs to be inverted, since it is literally covariance; 
  % the tikonov matrices below ("s") need not be inverted
  if invtype == 0
    rcov2 = inv(X2);
  elseif invtype == 1
    rcov2 = pinv(X2);
  elseif invtype == 2
    rcov2 = invillco(X2);
  elseif invtype == 3
    rcovF2 = factorize(X2);
    rcov2  = inverse(X2);
    rcov2  = rcov2 * eye(size(rcov2));
  elseif invtype == 4
    rcov2 = inverse_ridge_regression_matrix(X2,kmax);
  elseif invtype == 5
    rcov2 = inverse_minimum_eigenvalue_matrix_optim(X2,kmaxrange,sigminrange,'X2');
  end
  
  % Use following line for only Tikhonov reg THIS SHOULD NOT BE INVERTED
  l = get_l(driver.jacobian.numlays,1);    
  s = transpose(l)*l;
  
  %% now build the Tikhonov regularization block matrix, using "s"
  lenS = length(driver.jacobian.scalar_i);
  %% default : only column/stemp jacs, layer WV
  rc2 = blkdiag(zeros(lenS,lenS),driver.oem.alpha_water*s*reducefact);
  
  %% always assumes you want to fit for WV ... may want to keep T fixed 9so no fit) and also may not want to fit O3
  if isfield(driver.oem,'alpha_temp') & isfield(driver.oem,'alpha_ozone')
    rc2 = blkdiag(zeros(lenS,lenS),driver.oem.alpha_water*s*reducefact,driver.oem.alpha_temp*s*reducefact,driver.oem.alpha_ozone*s*reducefact);
  elseif isfield(driver.oem,'alpha_temp') & ~isfield(driver.oem,'alpha_ozone')
    rc2 = blkdiag(zeros(lenS,lenS),driver.oem.alpha_water*s*reducefact,driver.oem.alpha_temp*s*reducefact);
  elseif ~isfield(driver.oem,'alpha_temp') & isfield(driver.oem,'alpha_ozone')
    rc2 = blkdiag(zeros(lenS,lenS),driver.oem.alpha_water*s*reducefact,driver.oem.alpha_ozone*s*reducefact);
  end

  % if invtype == 3, rcov could be a class rather than double
  switch driver.oem.reg_type
    case 'reg_and_cov'
      r2 = rcov2 + rc2;
    case 'reg'
      r2 = rc2;
    case 'cov'
      r2 = rcov2;
    otherwise
      disp('Incorrect choice driver.oem.reg_type')
  end  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if length(intersect(iaSequential,210)) == 1
  reducefact = sqrt(10);
  reducefact = sqrt(10)/2;
  reducefact = sqrt(2);
  couplefact = 0.1;
  couplefact = 0.5;
  ioffsetWV = length(driver.jacobian.scalar_i) + length(driver.jacobian.water_i) - 3;
  ioffsetT  = length(driver.jacobian.scalar_i) + length(driver.jacobian.water_i) + length(driver.jacobian.temp_i) - 3;

  rcov210 = driver.oem.cov;
  X210 = driver.oem.cov;
  X210 = X210/(reducefact*reducefact); %% make unc^2 in parameters reducefact times 
  rcov210 = X210;

  boo = diag(rcov210);
  booT  = sqrt(boo(driver.jacobian.temp_i));
  booWV = sqrt(boo(driver.jacobian.water_i));
  for iiii = 1 : 4
    for jjjj = 1 : 4
      if iiii ~= jjjj
        rcov210(ioffsetWV+iiii,ioffsetT+jjjj) = couplefact*booT(iiii)*booWV(jjjj);
        rcov210(ioffsetT+iiii,ioffsetWV+jjjj) = couplefact*booT(iiii)*booWV(jjjj);
      end
    end
  end
  X210 = rcov210;

  % driver.oem.cov needs to be inverted, since it is literally covariance; 
  % the tikonov matrices below ("s") need not be inverted
  if invtype == 0
    rcov210 = inv(X210);
  elseif invtype == 1
    rcov210 = pinv(X210);
  elseif invtype == 2
    rcov210 = invillco(X210);
  elseif invtype == 3
    rcovF210 = factorize(X210);
    rcov210  = inverse(X210);
    rcov210  = rcov210 * eye(size(rcov210));
  elseif invtype == 4
    rcov210 = inverse_ridge_regression_matrix(X210,kmax);
  elseif invtype == 5
    rcov210 = inverse_minimum_eigenvalue_matrix_optim(X210,kmaxrange,sigminrange,'X210');
  end
  
  % Use following line for only Tikhonov reg THIS SHOULD NOT BE INVERTED
  l = get_l(driver.jacobian.numlays,1);    
  s = transpose(l)*l;
  
  %% now build the Tikhonov regularization block matrix, using "s"
  lenS = length(driver.jacobian.scalar_i);
  %% default : only column/stemp jacs, layer WV
  rc210 = blkdiag(zeros(lenS,lenS),driver.oem.alpha_water*s*reducefact);
  
  %% always assumes you want to fit for WV ... may want to keep T fixed 9so no fit) and also may not want to fit O3
  if isfield(driver.oem,'alpha_temp') & isfield(driver.oem,'alpha_ozone')
    rc210 = blkdiag(zeros(lenS,lenS),driver.oem.alpha_water*s*reducefact,driver.oem.alpha_temp*s*reducefact,driver.oem.alpha_ozone*s*reducefact);
  elseif isfield(driver.oem,'alpha_temp') & ~isfield(driver.oem,'alpha_ozone')
    rc210 = blkdiag(zeros(lenS,lenS),driver.oem.alpha_water*s*reducefact,driver.oem.alpha_temp*s*reducefact);
  elseif ~isfield(driver.oem,'alpha_temp') & isfield(driver.oem,'alpha_ozone')
    rc210 = blkdiag(zeros(lenS,lenS),driver.oem.alpha_water*s*reducefact,driver.oem.alpha_ozone*s*reducefact);
  end

  % if invtype == 3, rcov could be a class rather than double
  switch driver.oem.reg_type
    case 'reg_and_cov'
      r210 = rcov210 + rc210;
    case 'reg'
      r210 = rc210;
    case 'cov'
      r210 = rcov210;
    otherwise
      disp('Incorrect choice driver.oem.reg_type')
  end  

  if length(intersect(iaSequential,214)) == 1
    r214    = r210;
    rcov214 = rcov210;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if length(intersect(iaSequential,214)) == 1
  reducefact = sqrt(10);
  reducefact = sqrt(10)/2;
  reducefact = sqrt(2);
  couplefact = 0.1;
  couplefact = 0.5;
  ioffsetWV = length(driver.jacobian.scalar_i) + length(driver.jacobian.water_i) - 3;
  ioffsetT  = length(driver.jacobian.scalar_i) + length(driver.jacobian.water_i) + length(driver.jacobian.temp_i) - 3;

  rcov210 = driver.oem.cov;
  X210 = driver.oem.cov;
  X210 = X210/(reducefact*reducefact); %% make unc^2 in parameters reducefact times 
  rcov210 = X210;

  boo = diag(rcov210);
  booT  = sqrt(boo(driver.jacobian.temp_i));
  booWV = sqrt(boo(driver.jacobian.water_i));
  for iiii = 1 : 4
    for jjjj = 1 : 4
      if iiii ~= jjjj
        rcov210(ioffsetWV+iiii,ioffsetT+jjjj) = couplefact*booT(iiii)*booWV(jjjj);
        rcov210(ioffsetT+iiii,ioffsetWV+jjjj) = couplefact*booT(iiii)*booWV(jjjj);
      end
    end
  end
  X210 = rcov210;

  % driver.oem.cov needs to be inverted, since it is literally covariance; 
  % the tikonov matrices below ("s") need not be inverted
  if invtype == 0
    rcov210 = inv(X210);
  elseif invtype == 1
    rcov210 = pinv(X210);
  elseif invtype == 2
    rcov210 = invillco(X210);
  elseif invtype == 3
    rcovF210 = factorize(X210);
    rcov210  = inverse(X210);
    rcov210  = rcov210 * eye(size(rcov210));
  elseif invtype == 4
    rcov210 = inverse_ridge_regression_matrix(X210,kmax);
  elseif invtype == 5
    rcov210 = inverse_minimum_eigenvalue_matrix_optim(X210,kmaxrange,sigminrange,'X210');
  end
  
  % Use following line for only Tikhonov reg THIS SHOULD NOT BE INVERTED
  l = get_l(driver.jacobian.numlays,1);    
  s = transpose(l)*l;
  
  %% now build the Tikhonov regularization block matrix, using "s"
  lenS = length(driver.jacobian.scalar_i);
  %% default : only column/stemp jacs, layer WV
  rc210 = blkdiag(zeros(lenS,lenS),driver.oem.alpha_water*s*reducefact);
  
  %% always assumes you want to fit for WV ... may want to keep T fixed 9so no fit) and also may not want to fit O3
  if isfield(driver.oem,'alpha_temp') & isfield(driver.oem,'alpha_ozone')
    rc210 = blkdiag(zeros(lenS,lenS),driver.oem.alpha_water*s*reducefact,driver.oem.alpha_temp*s*reducefact,driver.oem.alpha_ozone*s*reducefact);
  elseif isfield(driver.oem,'alpha_temp') & ~isfield(driver.oem,'alpha_ozone')
    rc210 = blkdiag(zeros(lenS,lenS),driver.oem.alpha_water*s*reducefact,driver.oem.alpha_temp*s*reducefact);
  elseif ~isfield(driver.oem,'alpha_temp') & isfield(driver.oem,'alpha_ozone')
    rc210 = blkdiag(zeros(lenS,lenS),driver.oem.alpha_water*s*reducefact,driver.oem.alpha_ozone*s*reducefact);
  end

  % if invtype == 3, rcov could be a class rather than double
  switch driver.oem.reg_type
    case 'reg_and_cov'
      r210 = rcov210 + rc210;
    case 'reg'
      r210 = rc210;
    case 'cov'
      r210 = rcov210;
    otherwise
      disp('Incorrect choice driver.oem.reg_type')
  end  

  if length(intersect(iaSequential,214)) == 1
    r214    = r210;
    rcov214 = rcov210;
  end
end

