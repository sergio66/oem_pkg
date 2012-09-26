function driver = get_rates(driver)
% Read the rate file and pull in obs, cal or bias rates.

if length(strfind(driver.rateset.datafile,'iasiB1')) == 0
  %% looks like all spectral rates are in ONE file
  load(driver.rateset.datafile)
else
  bwah1 = load(driver.rateset.datafile);
  xdatafile = driver.rateset.datafile;
  ooh = strfind(driver.rateset.datafile,'iasiB1');
  xdatafile(ooh+5:ooh+5) = '2'
  bwah2 = load(xdatafile);

  b_obs = cat(2,bwah1.b_obs,bwah2.b_obs);
  b_cal = cat(2,bwah1.b_cal,bwah2.b_cal);
  b_bias = cat(2,bwah1.b_bias,bwah2.b_bias);

  b_err_obs = cat(2,bwah1.b_err_obs,bwah2.b_err_obs);
  b_err_cal = cat(2,bwah1.b_err_cal,bwah2.b_err_cal);
  b_err_bias = cat(2,bwah1.b_err_bias,bwah2.b_err_bias);

  clear bwah1 bwah2
  
end

% Either obs, bias, cal for the selected latitude
iibin = driver.iibin;

switch driver.rateset.ocb_set
  case 'bias'
     driver.rateset.rates = real(squeeze(b_bias(iibin,:,2))');
     driver.rateset.unc_rates = real(squeeze(b_err_bias(iibin,:,2))');
  case 'cal'
     driver.rateset.rates = real(squeeze(b_cal(iibin,:,2))');
     driver.rateset.unc_rates = real(squeeze(b_err_cal(iibin,:,2))');
  case {'obs','tracegas'}
     driver.rateset.rates = real(squeeze(b_obs(iibin,:,2))');
     driver.rateset.unc_rates = real(squeeze(b_err_obs(iibin,:,2))');
end

bad = find(driver.rateset.rates > driver.rateset.max);
if length(bad) > 0
  disp('resetting some bad input (max)');
  driver.rateset.rates(bad) = driver.rateset.max;
end
bad = find(driver.rateset.rates < driver.rateset.min);
if length(bad) > 0
  disp('resetting some bad input (min)');
  driver.rateset.rates(bad) = driver.rateset.min;
end



