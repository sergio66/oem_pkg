iTik = 1; %% order of tikonov = 0,1,2  DEFAULT

if length(iaSequential) == 1
  error('this is sequentail retrieval, so expect length(iaSequential) > 1')
end
if length(iaSequential) > 1 & iaSequential(1) == -1
  error('this is sequentail retrieval, so expect length(iaSequential) > 1 and you DO NOT EXPECT iaSequential(1) == -1')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if iSequential == -1
  %% remember, you expect this AFTER iSequential == 150,60,100 have been done, so don't change retrieval too much after that
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
  l = get_l(driver.jacobian.numlays,iTik);    
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

  r = r2;
  rcov = rcov2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif iSequential == 150 | iSequential == 60 | iSequential == 100
  X2 = driver.oem.cov;
  X2 = X2(iUseRetrParam,iUseRetrParam);
  reducefact = 1;
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
  l = get_l(driver.jacobian.numlays,iTik);    
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
  rc2 = rc2(iUseRetrParam,iUseRetrParam);

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

  r = r2;
  rcov = rcov2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif iSequential == 214
  reducefact = sqrt(10);
  reducefact = sqrt(10)/2;
  reducefact = sqrt(2);
  reducefact = 1;
  couplefact = 0.1;
  couplefact = 0.5;

  ioffsetST = 6;
  ioffsetWV = length(driver.jacobian.scalar_i) + length(driver.jacobian.water_i) - (iNXYZLay-1);
  ioffsetT  = length(driver.jacobian.scalar_i) + length(driver.jacobian.water_i) + length(driver.jacobian.temp_i) - (iNXYZLay-1);

  rcov210 = driver.oem.cov;
  X210 = driver.oem.cov;
  X210 = X210/(reducefact*reducefact); %% make unc^2 in parameters reducefact times 
  rcov210 = X210;

  boo = diag(rcov210);
  booST = sqrt(boo(driver.jacobian.scalar_i(6)));
  booT  = sqrt(boo(driver.jacobian.temp_i));
  booWV = sqrt(boo(driver.jacobian.water_i));
  for iiii = 1 : iNXYZLay
    for jjjj = 1 : iNXYZLay
      if iiii ~= jjjj
        rcov210(ioffsetWV+iiii,ioffsetT+jjjj) = couplefact*booT(iiii)*booWV(jjjj) * exp(-(abs(iiii-jjjj)));
        rcov210(ioffsetT+jjjj,ioffsetWV+iiii) = couplefact*booT(jjjj)*booWV(iiii) * exp(-(abs(iiii-jjjj)));

        rcov210(ioffsetST,ioffsetT+jjjj) = couplefact*booST*booT(jjjj) * exp(-(abs(0-jjjj)));
        rcov210(ioffsetT+jjjj,ioffsetST) = couplefact*booST*booT(jjjj) * exp(-(abs(0-jjjj)));

        rcov210(ioffsetST,ioffsetWV+iiii) = couplefact*booST*booWV(iiii) * exp(-(abs(0-iiii)));
        rcov210(ioffsetWV+iiii,ioffsetST) = couplefact*booST*booWV(iiii) * exp(-(abs(0-iiii)));
      end
    end
  end

  X210 = rcov210(iUseRetrParam,iUseRetrParam);

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
  l = get_l(driver.jacobian.numlays,iTik);    
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
  rc210 = rc210(iUseRetrParam,iUseRetrParam);

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

  r    = r210;
  rcov = rcov210;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
