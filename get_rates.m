function driver = get_rates(driver)
% Read the rate file and pull in obs, cal or bias rates.

load(driver.rateset.datafile)

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



