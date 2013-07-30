function driver = set_struct();

%---------------------------------------------------------------------------
% Set code debugger   
driver.debug = false;

%---------------------------------------------------------------------------
% Perform OEM fit
% Generallly true, but if you want to do LLS only, set this to false
driver.oem.dofit = true;

%---------------------------------------------------------------------------
% default is set up for AIRS rates, bin 20 
%    qst = 6 = [CO2 O3 N2O CH4 CFC stemp]
%    97 lays for Q1/T   
%---------------------------------------------------------------------------

% Which latitude bin
driver.iibin = 20;

% Debug dir
driver.debug_dir = '../Debug';

% Output filename
driver.filename = ['test1_' int2str(driver.iibin)];

% Observed rate file
driver.rateset.datafile = 'fitout_8year_v4_robust.mat';

% Covariance file for retrieved params Q WV T
driver.oem.cov_filename = 'cov_lls.mat';

% Covariance file for spectra : assume originally it DoesNotExist
driver.oem.spectralcov_filename = 'dne';

% sarta Forward Model error at all channels
driver.oem.sarta_error = 0.0;

% Apriori file
driver.oem.apriori_filename = 'apriori_zero';

% Lag-1 correlation file
driver.rateset.ncfile   = 'all_lagcor.mat';
driver.rateset.ncfile   = 'dne';

% Jacobian file: f = 2378x1 and M_TS_jac_all = 36x2378x200
driver.jacobian.filename = 'M_TS_jac_all.mat';

% Select channels to fit
load goodchan.mat
driver.jacobian.chanset = goodchan;   % 1500 channels

% do we get rid of detected spikes ( <0 is NO, >0 is yes)
driver.rateset.despike = -1;

%---------------------------------------------------------------------------
% Do one of [obs][cal][biases] 
driver.rateset.ocb_set  = 'obs';

% when you read in the delta(BTs) you can adjust them by multiplying by a scale factor
% this affect the Se matrix, and therefore the reported uncertainties in retrieved params
driver.oem.adjust_spectral_errorbars = 1;

% Remove known rates?
driver.rateset.adjust  = false;
driver.rateset.adjust_index  = [1   103 200];
driver.rateset.adjust_values = [2.2 0.01 0.01];  %% physical units

driver.jacobian.numQlays   = 1;    %% have N lays of temp jacs; how many gases we want to retrieve profiles
                                   %% must be at least 1 (which we assume is water, though not necessarily!)
driver.jacobian.gasID_Qlays = [1]; %% ID of gases whose profile you are retrieving [1] is minimum necessary
driver.jacobian.numlays    = 97;

% Choose QST variables to fit
% note preference : CO2 first, SST last ... others in between!
driver.jacobian.qstnames = {'CO2' 'O3' 'N2O' 'CH4' 'CFC11' 'stemp'};
driver.jacobian.qstYesOrNo = [1    1     1     1     1       1];

% Choose wv and t have maximum 97 layers (1=toa,97 = gnd)
driver.jacobian.Q1jacindex  = 1:97;
driver.jacobian.tjacindex   = 1:97;

% lambda used if just multiplying diagonal
driver.oem.lambda = 2;

% This lambda multiplying diag only (true)? or entire matrix (false)? 
driver.oem.diag_only = false;

% SARTA forward model and other "representation" errors
driver.sarta_error = 0.0; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% how many times are we looping
driver.oem.nloop = 1;

% are we doing rates (+1) or regular spectra (-1)
driver.oem.rates = +1;
% these are needed if doing regular spectra ie if driver.oem.rates = -1;
  driver.oem.klayers = '/asl/packages/klayersV205/BinV201/klayers_airs';
  driver.oem.sarta   = '/asl/packages/sartaV108_PGEv6/Bin/sarta_airs_PGEv6_postNov2003';
  driver.oem.sarta   = '/asl/packages/sartaV108/Bin/sarta_apr08_m140_iceaggr_waterdrop_desertdust_slabcloud_hg3_wcon_nte';
  driver.oem.headstruct  = [];      %% this should be the starting head structure
  driver.oem.profstruct  = [];      %% this should be the starting profile structure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NOTE THE DEFAULTS ARE FOR AIRS

load airs_test_chanset.mat
driver.jacobian.chanset = chanset;
clear chanset

%% these are from strow_stmNov2011_dobs20.mat, where Larrabee says he gets
%%    good CH4 rates
driver.jacobian.filename = 'M_TS_jac_all.mat';

driver.rateset.datafile = 'fitout_8year_v4_robust.mat';
driver.rateset.ocb_set  = 'obs';
driver.rateset.adjust   = 0;

>>>>>>
%% these are the geophysical covariances and/or smoothing params
driver.oem.regularizationVScovariances = 'R';
  %% this means set the SMOOTHING L1/L2 params
  driver.oem.cov_filename     = 'cov_lls.mat';
  driver.oem.apriori_filename = 'apriori_zero';
  driver.oem.diag_only        = 0;
  driver.oem.lambda           = 1;
  driver.oem.lambda_qst       = 0.1;
  driver.oem.lambda_Q1        = 10;
  driver.oem.lambda_temp      = 1;

%% we've settled on geophysical covariance for clear sky rates
driver.oem.regularizationVScovariances = 'C';
  %% but if smoothVSregularization = 'c'
  %% note these are PHYSICAL units eg sigma_temp_stratVALUE is in KELVIN (or KELVIN/YR)
  %% so eg Andy Tangborn finds after normalization driver.oem.sigma.temp_strat_VALUE  = 4;
  %% is a "good value" which means we default set it to 4*0.01 (0.01 is qrenorm for T(z))
  %% strat temp
    driver.oem.sigma.temp_strat_VALUE  = 4.0*0.01;       %% sigsqr
    driver.oem.sigma.temp_strat_TOPLAY = 01;  %% start layer
    driver.oem.sigma.temp_strat_BOTLAY = 49;  %% stop layer
  %% trop temp
    driver.oem.sigma.temp_trop_VALUE  = 0.5*0.01;       %% sigsqr
    driver.oem.sigma.temp_trop_TOPLAY = 50;  %% start layer
    driver.oem.sigma.temp_trop_BOTLAY = 97;  %% stop layer
    driver.oem.sigma.temp_trop_TOPLAY = driver.oem.sigma.temp_strat_BOTLAY + 1; %% start layer
    driver.oem.sigma.temp_trop_BOTLAY = driver.jacobian.numlays;                %% stop layer
  %% strat GAS 1
    driver.oem.sigma.hum_strat_VALUE  = 1.5*0.01;       %% sigsqr
    driver.oem.sigma.hum_strat_TOPLAY = 01;  %% start layer
    driver.oem.sigma.hum_strat_BOTLAY = 49;  %% stop layer
  %% trop GAS 1
    driver.oem.sigma.hum_trop_VALUE  = 1.0*0.01;       %% sigsqr
    driver.oem.sigma.hum_trop_TOPLAY = 50;  %% start layer
    driver.oem.sigma.hum_trop_BOTLAY = 97;  %% stop layer
    driver.oem.sigma.hum_trop_TOPLAY = driver.oem.sigma.hum_strat_BOTLAY + 1; %% start layer
    driver.oem.sigma.hum_trop_BOTLAY = driver.jacobian.numlays;                %% stop layer
  %% length_correlation for T(z) and WV(z) .. note this is in terms of INDEX
    driver.oem.sigma.l_c = 2.4;
  %% QST 1 ..6 values
    driver.oem.sigma.qst(1) = 4*2.20;  %% co2
    driver.oem.sigma.qst(2) = 1*0.01;  %% o3
    driver.oem.sigma.qst(3) = 4*1.00;  %% n2o
    driver.oem.sigma.qst(4) = 1*5.00;  %% ch4 
    driver.oem.sigma.qst(5) = 1*1.00;  %% cfc11
    driver.oem.sigma.qst(6) = 1*0.10;  %% Stemp

