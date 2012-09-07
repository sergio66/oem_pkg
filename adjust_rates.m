function driver = adjust_rates(driver,m_ts_jac);

% Removes "known" components of the rate spectra

coeffs = driver.rateset.adjust_values ./ driver.qrenorm(driver.rateset.adjust_index);

known_portion = zeros(size(driver.rateset.rates));
for ii = 1 : length(coeffs)
  ix = driver.rateset.adjust_index(ii);
  known_portion = known_portion + coeffs(ii)*m_ts_jac(:,ix);
end
driver.rateset.rates = driver.rateset.rates - known_portion;
