function driver = set_struct();

%---------------------------------------------------------------------------
% Set code debugger   
driver.debug = false;

%---------------------------------------------------------------------------
% Perform OEM fit
% Generallly true, but if you want to do LLS only, set this to false
driver.oem.dofit = true;

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

% Remove known rates?
driver.rateset.adjust  = false;
driver.rateset.adjust_index  = [1   103 200];
driver.rateset.adjust_values = [2.2 0.01 0.01];  %% physical units

% Choose variables to fit
driver.jacobian.co2   = true;
driver.jacobian.o3    = true;
driver.jacobian.n2o   = true;
driver.jacobian.ch4   = true;
driver.jacobian.cfc11 = true;
driver.jacobian.stemp = true;
% wv and t have maximum 97 layers (1=toa,97 = gnd)
driver.jacobian.wvjacindex  = 1:97;
driver.jacobian.tjacindex   = 1:97;

% lambda used if just multiplying diagonal
driver.oem.lambda = 2;

% This lambda multiplying diag only (true)? or entire matrix (false)? 
driver.oem.diag_only = false;

% SARTA forward model and other "representation" errors
driver.sarta_error = 0.0; 

