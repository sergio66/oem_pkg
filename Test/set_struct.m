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

% Covariance file
driver.oem.cov_filename = 'cov_lls.mat';

% Apriori file
driver.oem.apriori_filename = 'apriori_zero';

% Lag-1 correlation file
driver.rateset.ncfile   = 'all_lagcor.mat';

% Jacobian file: f = 2378x1 and M_TS_jac_all = 36x2378x200
driver.jacobian.filename = 'M_TS_jac_all.mat';

% Select channels to fit
load goodchan.mat
driver.jacobian.chanset = goodchan;   % 1500 channels

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

driver.jacobian.numQlays   = 1;   %% have N lays of temp jacs; how many gases we want to retrieve profiles
                                  %% must be at least 1 (which we assume is water, though not necessarily!)
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

% are we having another gas other than wv (eg HDO???)
driver.oem.othergases = -1;   %% -1 means NO, +X means yes, X more gases

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NOTE THE DEFAULTS ARE FOR AIRS

load airs_test_chanset.mat
driver.jacobian.chanset = chanset;
clear chanset

%% these are from strow_stmNov2011_dobs20.mat, where Larrabee says he gets
%%    good CH4 rates
driver.jacobian.filename = 'M_TS_jac_all.mat';

driver.rateset.datafile = 'fitout_8year_v4_robust.mat';
driver.rateset.ncfile   = 'all_lagcor.mat';
driver.rateset.ocb_set  = 'obs';
driver.rateset.adjust   = 0;

driver.oem.cov_filename     = 'cov_lls.mat';
driver.oem.apriori_filename = 'apriori_zero';
driver.oem.diag_only        = 0;
driver.oem.lambda           = 1;
%% >>>> orig
%% driver.oem.lambda_gas       = 0.1;
%% driver.oem.lambda_water     = 10;
%% driver.oem.lambda_temp      = 1;
%%% >>> new
driver.oem.lambda_qst       = 0.1;
driver.oem.lambda_Q1        = 10;
driver.oem.lambda_temp      = 1;

