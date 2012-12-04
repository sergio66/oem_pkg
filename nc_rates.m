function ncerrors =  nc_rates(driver);

% Correct rate uncertainties with 1-lag autocorrelation from ncfile

junkname = driver.jacobian.filename;

nc = ones(1,2378);
if findstr(junkname,'IASI')
  nc = ones(1,8642);
end

%% this loop is mainly for AIRS; still need to get nc for IASI
if ~strcmp('dne',driver.rateset.ncfile)
  fprintf(1,'   using nc errors from %s \n',driver.rateset.ncfile)
  load(driver.rateset.ncfile)
  iibin = driver.iibin;
  switch driver.rateset.ocb_set
    case 'bias'
      nc = (1+lagcor_bias_anom(iibin,:))./(1-lagcor_bias_anom(iibin,:));
      nc = sqrt(nc);
    case 'cal'
      nc = (1+lagcor_cal_anom(iibin,:))./(1-lagcor_cal_anom(iibin,:));
      nc = sqrt(nc);
    case {'obs','tracegas'}
      nc = (1+lagcor_obs_anom(iibin,:))./(1-lagcor_obs_anom(iibin,:));
      nc = sqrt(nc);
  end
end

if findstr(junkname,'IASI')
  %% this is for IASI
  junk = driver.rateset.unc_rates;
  ncerrors = junk;
  disp('nc errors not yet coded up for IASI, default to unc_rates')
else
  %% this is for AIRS
  ncerrors = real(nc' .* driver.rateset.unc_rates);   
end
