function r = regularization_multiplier(r,driver);

%% this tweaks the 200x 200 parameter covariance matrix
%% according to the various "lambda"
%% called by rodgers.m

if driver.oem.diag_only
  % Paul's suggestion to access the diagonal
  r(1:size(r,1)+1:end) = r(1:size(r,1)+1:end) * driver.oem.lambda;

else
  %%%%%%%%%%%%%%%%%%%% this is for QST ie the SINGLE COLUMN gases and SST
  % Multiply entire matrix, but in blocks of column(QST), qlays(Q1 .. QN), templays
  if length(driver.oem.lambda_qst) == 1
    %% multiply all R(i,j) of the qst part of the matrix, by one number
    r(driver.jacobian.iqst,driver.jacobian.iqst)     = r(driver.jacobian.iqst,driver.jacobian.iqst) * driver.oem.lambda_qst;
  elseif length(driver.oem.lambda_qst) == length(driver.jacobian.iqst)
    %% multiply all R(i,i) of the qst part of the matrix, ,* specified diagnol, making this diagnol
    r(driver.jacobian.iqst,driver.jacobian.iqst)     = r(driver.jacobian.iqst,driver.jacobian.iqst) .* diag(driver.oem.lambda_qst);
  elseif length(driver.oem.lambda_qst) > length(driver.jacobian.iqst)
    %% matrix multiply all R(i,j) of the qst part of the matrix, by specified diagnol
    junk = driver.oem.lambda_qst(1:length(driver.jacobian.iqst));
    junk = diag(junk);
    r(driver.jacobian.iqst,driver.jacobian.iqst)     = r(driver.jacobian.iqst,driver.jacobian.iqst) * junk;
  end

  %%%%%%%%%%%%%%%%%%%% this is for Q(z) ie the MULTICOLUMN gases
  for ii = 1 : driver.jacobian.numQlays
    junk = ['lala = driver.oem.lambda_Q' num2str(ii) ';']; eval(junk);
    junk = ['mama = driver.jacobian.iQ' num2str(ii) ';'];  eval(junk);
    if length(lala) == 1
      %% multiply all R(i,j) of the iQ part of the matrix, by one number
      r(mama,mama) = r(mama,mama) * lala;
    elseif length(lala) == length(mama)
      %% multiply all R(i,i) of the iQ part of the matrix, .* specified diagnol, making this diagnol
      r(mama,mama) = r(mama,mama)  .* diag(lala);
    elseif length(lala) > length(mama)
      %% multiply all R(i,j) of the iQ part of the matrix, by specified diagnol
      junk = lala(1:length(mama));
      junk = diag(junk);
      r(mama,mama) = r(mama,mama)  * (junk);
    end
  end

  %%%%%%%%%%%%%%%%%%%% this is for T(z) ie the MULTICOLUMN temperature
  lala = driver.oem.lambda_temp;
  mama = driver.jacobian.itemp;
  if length(driver.oem.lambda_temp) == 1
    %% multiply all R(i,j) of the T part of the matrix, by one number
    r(driver.jacobian.itemp,driver.jacobian.itemp)   = r(driver.jacobian.itemp,driver.jacobian.itemp) * driver.oem.lambda_temp;
  elseif length(lala) == length(mama)
    %% multiply all R(i,i) of the T part of the matrix, .* specified diagnol, making this diagnol
    r(driver.jacobian.itemp,driver.jacobian.itemp)   = r(driver.jacobian.itemp,driver.jacobian.itemp) .* diag(driver.oem.lambda_temp);
  elseif length(lala) > length(mama)
    %% multiply all R(i,j) of the iQ part of the matrix, by specified diagnol
    junk = lala(1:length(mama));
    junk = diag(junk);
    r(mama,mama) = r(mama,mama)  * (junk);
  end

  % Now add diagonal, but only for T and water
  n = length(driver.jacobian.iQ1);
  m = length(driver.jacobian.itemp);
  ti = driver.jacobian.itemp;
  for i=1:m
    r(ti(i),ti(i)) = r(ti(i),ti(i)) + driver.oem.lambda;
  end

  for ii = 1 : driver.jacobian.numQlays
    junk = ['qi = driver.jacobian.iQ' num2str(ii) ';']; eval(junk);
    for i=1:n
      r(qi(i),qi(i)) = r(qi(i),qi(i)) + driver.oem.lambda;
    end
  end

  %% fine add it on for QST as well (Aug 2012)
  n = length(driver.jacobian.qstYesOrNo);
  for i=1:n
    r(i,i) = r(i,i) + driver.oem.lambda;
  end

end
