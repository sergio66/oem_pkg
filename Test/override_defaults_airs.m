function driver = override_defaults(driver,ix);

% can use the package eg as
% for ix=1:36                                  
%   ix
%   clear driver; run_retrieval; end

load airs_test_chanset.mat
driver.jacobian.chanset = chanset;

%% these are limits for good and bad input SPECTRA
driver.rateset.max = 320;
driver.rateset.min = 180;
%% these are limits for good and bad input SPECTRAL RATES
driver.rateset.max = +1;
driver.rateset.min = -1;

%% these are from strow_stmNov2011_dobs20.mat, where Larrabee says he gets
%%    good CH4 rates
driver.jacobian.filename = 'M_TS_jac_all.mat';

driver.rateset.datafile = 'fitout_8year_v4_robust.mat';
driver.rateset.ncfile   = 'all_lagcor.mat';
driver.rateset.ocb_set  = 'obs';
driver.rateset.adjust   = 0;

%%% larrabee used these in the file he gave me, but I get bad trace gas rates
%%% they give decent WV/T rates
driver.oem.cov_filename     = 'cov_lls.mat';
driver.oem.apriori_filename = 'apriori_zero';
driver.oem.diag_only        = 0;
driver.oem.lambda           = 1;
driver.oem.lambda_gas       = 0.1;
driver.oem.lambda_water     = 10;
driver.oem.lambda_temp      = 1;

%%% these give decent CO2, N2O rates!!!!! SAVE THIS!!!!
%%% these give decent CO2, N2O rates!!!!! SAVE THIS!!!!
%driver.oem.lambda           = 0.01;
%driver.oem.lambda_gas       = 0.01;
%driver.oem.lambda_water     = 0.01;
%driver.oem.lambda_temp      = 0.01;
%%% these give decent CO2, N2O rates!!!!! SAVE THIS!!!!
%%% these give decent CO2, N2O rates!!!!! SAVE THIS!!!!

%%%% not bad for all (trace gases and T,WV)
driver.oem.diag_only        = 0;
driver.oem.lambda           = 0.01;
driver.oem.lambda_gas       = 0.1;
driver.oem.lambda_water     = 10;
driver.oem.lambda_temp      = 10;
driver.oem.lambda_water     = 5;
driver.oem.lambda_temp      = 5;
%%%% not bad for all (trace gases and T,WV)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

driver.iibin = ix;
driver.filename = ['../Output/testx_' int2str(driver.iibin)];
