function driver = override_defaults(driver,ix);

% can use the package eg as
% for ix=1:36                                  
%   ix
%   clear driver; run_retrieval; end

%% used for LLS in Aug 2012, for Utah meeting

driver.oem.apriori_filename = 'apriori_zero';

addpath /home/sergio/MATLABCODE
y = instr_chans('iasi');
indB = find(y <= 2205);
driver.jacobian.chanset = indB;

%% these are limits for good and bad input SPECTRA
driver.rateset.max = 320;
driver.rateset.min = 180;
%% these are limits for good and bad input SPECTRAL RATES
driver.rateset.max = +1;
driver.rateset.min = -1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% orig 6_97_97 stuff
driver.jacobian.filename = ...
  '/home/sergio/MATLABCODE/TMP_RATES_Fit_pkg/Aux_jacs_IASI/AUG21_2011/all_kcarta_jacsB1.mat';
driver.jacobian.qstnames = {'CO2' 'O3' 'N2O' 'CH4' 'CFC11' 'stemp'};
driver.jacobian.qstYesOrNo = [1     1   1     1      1       1];
driver.jacobian.numQlays   = 1;   %% have N lays of temp jacs; how many gases we want to retrieve profiles
                                  %% must be at least 1 (for water)
driver.jacobian.numlays    = 97;

driver.rateset.datafile = '/home/sergio/MATLABCODE/TMP_RATES_Fit_pkg/Cluster_IASI/strow_fake36_Aug2012_sergio_iasi_rates.mat';
driver.jacobian.filename = '/home/sergio/MATLABCODE/TMP_RATES_Fit_pkg/Cluster_IASI/strow_fake36_Aug2012_IASI_sergio_jacB1.mat';

% example setup if you want to individualize the lambdas
% driver.oem.lambda_qst       = [0.1 0.2 0.3 0.4 0.5 0.6]*0.10*10;
% driver.oem.lambda_Q1        = [ones(1,10) ones(1,30)*2 ones(1,57)*3]*100;
% driver.oem.lambda_temp      = [ones(1,10)*1 ones(1,30)*5 ones(1,57)*1]*100;

% keep lambdas constant
driver.oem.lambda        = 0.01;
driver.oem.diag_only     = 0;
driver.oem.lambda_qst    = 0.1;
driver.oem.lambda_Q1     = 100;
driver.oem.lambda_temp   = 50;

driver.jacobian.chanset = indB;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
driver.iibin = ix;
driver.filename = ['../Output/testx_' int2str(driver.iibin)];

