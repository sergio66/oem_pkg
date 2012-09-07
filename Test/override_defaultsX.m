function driver = override_defaultsX(driver);

disp('in DEBUG mode : get rid of override_defaultsX in run_retrieval.m')
disp('in DEBUG mode : get rid of override_defaultsX in run_retrieval.m')
disp('in DEBUG mode : get rid of override_defaultsX in run_retrieval.m')

load /strowdata1/shared/sergio/MATLABCODE/TMP_RATES_Fit_pkg/Cluster/strow_stmNov2011_dobs20
driver.jacobian.chanset = dobs20.jacobian.chanset;

%% these are from strow_stmNov2011_dobs20.mat, where Larrabee says he gets
%%    good CH4 rates
driver.jacobian.filename = 'M_TS_jac_all.mat';

driver.rateset.datafile = 'fitout_8year_v4_robust.mat'
driver.rateset.ncfile   = 'all_lagcor.mat';
driver.rateset.ocb_set  = 'obs';
driver.rateset.adjust   = 0;

driver.oem.cov_filename     = 'cov_lls.mat'
driver.oem.apriori_filename = 'apriori_zero';
driver.oem.diag_only        = 0;
driver.oem.lambda           = 1;
driver.oem.lambda_gas       = 0.1;
driver.oem.lambda_water     = 10;
driver.oem.lambda_temp      = 1;

