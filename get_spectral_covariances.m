function e0 = get_spectral_covariances(driver,ncerrors,inds);

%% this function builds the 2378 x 2378 spectral covs 
%% if driver.oem.spectralcov_filename == 'dne' then do not read in files else
%%  read in cov matrix and renormalize, based on ncerrors
%% called by rodgers.m

% Obs errors set with serial correlation correction
xe0 = diag(ncerrors(inds));

if ~strcmp('dne',driver.oem.spectralcov_filename)
  %% see Test/set_struct.m as this sets default name to 'dne'
  spectral_cov = load(driver.oem.spectralcov_filename);
  dacov = spectral_cov.dacov;
  %% make sure there are some entries here!!!! ie not all zeros
  if nansum(diag(abs(dacov))) > 0
    fprintf(1,'   using spectral covs from %s \n',driver.oem.spectralcov_filename)
    e0 = dacov(inds,inds);

    %% now need to normalize
    xdiag = diag(xe0);
    xij = xdiag * xdiag';
    thediag = diag(e0);
    oo = find(isnan(e0(:))); e0(oo) = 0;
    oo = find(isinf(e0(:))); e0(oo) = 0;
    %whos e0 xe0 xdiag
    e0 = e0 * nanmean(xdiag)/nanmean(diag(e0));
  else
    e0 = xe0;
    fprintf(1,'  Warning : not using spectral covs from %s as they are NaNs \n',driver.oem.spectralcov_filename)
  end
else  %% just use ncerrors
  e0 = xe0;
end
