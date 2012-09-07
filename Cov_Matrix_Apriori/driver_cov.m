driver = get_rates_structure;

if driver.cov_choice == 0
  xcov = find_cov_simple_anom(driver);
elseif driver.cov_choice == 1
  xcov = find_cov_simple(driver);
elseif driver.cov_choice == 2
  xcov = find_cov_deriv(driver);
elseif driver.cov_choice == 3
  xcov = find_cov_damp(driver);
elseif driver.cov_choice == 4
  [xcov,ERArates,ERAstd,ERArenorm,cutoffs] = find_cov_damp2(driver,1);
elseif driver.cov_choice == 5
  xcov = find_cov_derivoperators(driver);
else
  error('need driver.cov_choice between 1 and 2')
end

cov           = xcov * driver.scale;
savecov.input = driver;

saver = ['save ' driver.covfilename ' savecov cov'];
eval(saver);