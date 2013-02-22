function e0 = get_spectral_covariances(driver,ncerrors,inds);

%% if driver.oem.spectralcov_filename == 'dne' then 
%%   use only ncerrors to make a diaganol 2378x2378 matrix
%% else 
%%   read in cov matrix and renormalize, based on ncerrors, to build
%%   the 2378 x 2378 spectral covs 
%% end

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

    % whos e0 xe0 xdiag
    % figure(1); clf; imagesc(e0); colorbar

    e0 = e0 * nanmean(xdiag)/nanmean(diag(e0));

    % figure(1); clf; plot(diag(e0));  axis([1 length(inds) 0 3e-3]); grid; title('cov(new)')
    % figure(2); clf; imagesc(e0);     colorbar; title('cov(new)')
    % figure(3); clf; imagesc(xe0);    caxis([0 1e-3]);               colorbar; title('cov(ncerrors)');
    % figure(4); clf; plot(diag(xe0)); axis([1 length(inds) 0 3e-3]); grid; title('cov(ncerrors)')    
    % pause

  else
    e0 = xe0;
    fprintf(1,'  Warning : not using spectral covs from %s as they are NaNs \n',driver.oem.spectralcov_filename)
  end
else  %% simple !!! just use diag(ncerrors)
  e0 = xe0;
end
