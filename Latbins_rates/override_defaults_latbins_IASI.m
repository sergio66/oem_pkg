function driver = override_defaults(driver,ix);

% can use the package eg as
% for ix=1:36                                  
%   ix
%   clear driver; run_retrieval; end

driver.oem.apriori_filename = 'apriori_zero';
addpath /home/sergio/MATLABCODE
y = instr_chans('iasi');
ind = find(y <= 2205);
ind = (1:4:length(ind));
ind = ind(ind > 200);

ind = find((y >= 725 & y <= 800) | (y >= 930 & y <= 980) | ...
           (y >= 1250 & y <= 1500) | (y >= 2100 & y <= 2250)); 
ind = ind(1:2:length(ind));
ind = find((y >= 725 & y <= 800) | (y >= 930 & y <= 980) | ...
           (y >= 1250 & y <= 1650) | (y >= 2100 & y <= 2250)); 
ind = ind(1:2:length(ind)); indA = ind;

driver.jacobian.chanset = (1:4230);
driver.jacobian.chanset = (1:4:8200);
dfs = '/strowdata1/shared/sergio/MATLABCODE/DFS_OPTIMUMCHANS_RODGERS/';
c1 = load([dfs 'iasi_optchannelZ.mat']); %% o3
c2 = load([dfs 'iasi_optchannelM.mat']); %% ch4
c3 = load([dfs 'iasi_optchannelW.mat']); %% h2o
c4 = load([dfs 'iasi_optchannelT.mat']); %% t
c5 = unique([c1.iaFreqX; c2.iaFreqX; c3.iaFreqX; c4.iaFreqX]);
driver.jacobian.chanset = c5(c5 > 422);  %% so we start at 750 cm-1
%driver.jacobian.chanset = c5(c5 > 262);  %% so we start at 710 cm-1
%blonk = find(y < 980 | y > 1100);
%driver.jacobian.chanset = intersect(driver.jacobian.chanset,blonk);
driver.jacobian.chanset = [(300:5:4230)  582:592];
driver.jacobian.chanset = [[282:592] [1342:1820] [2642:3040]];
driver.jacobian.chanset = [[282:840] [1342:1820] [2662:3240]];
driver.jacobian.chanset = [[282:840] [1342:1820] [2602:3240]];
driver.jacobian.chanset = [[282:780] [1342:1740] [2602:3240]];
%driver.jacobian.chanset = c4.iaFreqX;
%  driver.jacobian.chanset = c4.iaFreqX(c4.iaFreqX > 422);  %%start at 750 cm-1
%  driver.jacobian.chanset = c4.iaFreqX(c4.iaFreqX > 262);  %%start at 710 cm-1
driver.jacobian.chanset = sort(unique(driver.jacobian.chanset));

indB = find(y <= 2205);
indC = find(y <= 2205 | y >= 2405);

driver.jacobian.chanset = indA;
driver.jacobian.chanset = indB;
driver.jacobian.chanset = indC;

xyz = '/strowdata1/shared/sergio/MATLABCODE/RATES_TARO/MAT/';
%% xyz = [xyz 'overocean_gsx_1dayV1_ecmwf2_iasi_lays_spanday01_avgL1Brates_robust_Aug27_2011.mat'];
xyz = [xyz 'overocean_gsx_1dayV1_ecmwf2_iasi_lays_spanday01_avgL1Brates_robust_Aug31_2011.mat'];
%% xyz = [xyz 'overocean_gsx_1dayV1_ecmwf2_iasi_lays_spanday01_avgL1Brates_robust_Feb1_2012.mat'];
driver.rateset.datafile = xyz;

%% by default, this should be for B1; B2 will automatically be read in
driver.jacobian.filename = ...
  '/home/sergio/MATLABCODE/TMP_RATES_Fit_pkg/Aux_jacs_IASI/AUG21_2011/all_kcarta_jacsB1.mat';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%driver.rateset.datafile =  '/home/sergio/MATLABCODE/TMP_RATES_Fit_pkg/Cluster_IASI/avg_strow_fake36_Aug2012_sergio_iasi_rates.mat';
driver.rateset.datafile = '/home/sergio/MATLABCODE/TMP_RATES_Fit_pkg/Cluster_IASI/strow_fake36_Aug2012_sergio_iasi_rates.mat';
driver.jacobian.filename = '/home/sergio/MATLABCODE/TMP_RATES_Fit_pkg/Cluster_IASI/strow_fake36_Aug2012_IASI_sergio_jacB1.mat';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% larrabee used these in the file he gave me, but I get bad trace gas rates
driver.oem.cov_filename     = 'cov_lls.mat'
driver.oem.apriori_filename = 'apriori_zero';
driver.oem.diag_only        = 0;
driver.oem.lambda           = 1;
driver.oem.lambda_qst       = 0.1;
driver.oem.lambda_Q1        = 10;
driver.oem.lambda_temp      = 1;

%%% these give decent CO2, N2O rates!!!!! (AIRS)
driver.oem.lambda           = 0.01 * 1;
driver.oem.lambda_qst       = 0.01 * 1;
driver.oem.lambda_Q1        = 0.01 * 1;
driver.oem.lambda_temp      = 0.01 * 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
driver.oem.apriori_filename = 'apriori_one';
driver.oem.diag_only        = 0;
driver.oem.lambda           = 0.01;
driver.oem.lambda_qst       = 0.1;
driver.oem.lambda_Q1        = 50;
driver.oem.lambda_temp      = 50;
driver.oem.lambda_Q1        = 0.1;
driver.oem.lambda_temp      = 0.1;

%% used for LLS in Aug 2012, for Utah meeting
driver.oem.lambda_Q1        = 10;
driver.oem.lambda_temp      = 1;
driver.oem.lambda_Q1        = 50;
driver.oem.lambda_temp      = 10;
driver.oem.lambda_Q1        = 100;
driver.oem.lambda_temp      = 50;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

driver.oem.lambda_qst       = [0.1 0.2 0.3 0.4 0.5 0.6]*0.10*10;
driver.oem.lambda_Q1        = [ones(1,10) ones(1,30)*2 ones(1,57)*3]*100;
driver.oem.lambda_temp      = [ones(1,10)*1 ones(1,30)*5 ones(1,57)*1]*100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% orig 6_97_97 stuff
driver.jacobian.filename = ...
  '/home/sergio/MATLABCODE/TMP_RATES_Fit_pkg/Aux_jacs_IASI/AUG21_2011/all_kcarta_jacsB1.mat';
driver.jacobian.qstnames = {'CO2' 'O3' 'N2O' 'CH4' 'CFC11' 'Stemp'};
driver.jacobian.qstYesOrNo = [1     1   1     1      1       1];
driver.jacobian.numQlays   = 1;   %% have N lays of temp jacs; how many gases we want to retrieve profiles
                                  %% must be at least 1 (for water)
driver.jacobian.numlays    = 97;
driver.oem.lambda_qst       = [0.1 0.2 0.3 0.4 0.5 0.6]*0.10*10;
driver.oem.lambda_Q1        = [ones(1,10) ones(1,30)*2 ones(1,57)*3]*100;
driver.oem.lambda_temp      = [ones(1,10)*1 ones(1,30)*5 ones(1,57)*1]*100;

driver.rateset.datafile =  '/home/sergio/MATLABCODE/TMP_RATES_Fit_pkg/Cluster_IASI/avg_strow_fake36_Aug2012_sergio_iasi_rates.mat';
driver.jacobian.filename = '/home/sergio/MATLABCODE/TMP_RATES_Fit_pkg_General/Cluster_IASI_General/avg_strow_fake36_Aug2012_IASI_sergio_jacfakeB1.mat';
driver.rateset.datafile = '/home/sergio/MATLABCODE/TMP_RATES_Fit_pkg/Cluster_IASI/strow_fake36_Aug2012_sergio_iasi_rates.mat';
driver.jacobian.filename = '/home/sergio/MATLABCODE/TMP_RATES_Fit_pkg/Cluster_IASI/strow_fake36_Aug2012_IASI_sergio_jacB1.mat';

driver.oem.diag_only        = 0;
driver.oem.lambda           = 0.01;
driver.oem.lambda_qst       = 0.1;
driver.oem.lambda_Q1     = 0.1;
driver.oem.lambda_temp      = 0.1;

driver.oem.lambda_Q1     = 10;
driver.oem.lambda_temp      = 1;
driver.oem.lambda_Q1     = 50;
driver.oem.lambda_temp      = 10;
driver.oem.lambda_Q1     = 100;
driver.oem.lambda_temp      = 50;

driver.jacobian.chanset = indB;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{

%% new 10_97_97 stuff
driver.jacobian.filename = '../Aux_jacs_IASI_General/AUG30_2012/all_kcarta_jacsB1_10_97_97.mat';
driver.jacobian.filename = '../Aux_jacs_IASI_General/AUG30_2012/all_kcarta_jacsB1_10_97_97.mat';
driver.jacobian.filename = '../Cluster_IASI_General/strow_fake36_Aug2012_IASI_sergio_jac_10_97_97fakeB1.mat';
driver.jacobian.filename = '../Cluster_IASI_General/avg_strow_fake36_Aug2012_IASI_sergio_jac_10_97_97fakeB1.mat';
driver.jacobian.qstnames = {'CO2trop' 'CO2strat' 'O3trop' 'O3strat' 'N2O' 'CO' 'CH4' 'CFC11' 'HDO' 'stemp'};
driver.jacobian.qstYesOrNo = [1       1          1        1         1     1     1    1       1     1];
driver.jacobian.numQlays   = 1;   %% have N lays of temp jacs; how many gases we want to retrieve profiles
                                  %% must be at least 1 (for water)
driver.jacobian.numlays    = 97;
%% just need a 204x204 for this
driver.oem.apriori_filename = '../Aux_jacs_AIRS_General/AUG30_2012/apriori204.mat';
driver.oem.cov_filename     = '../Aux_jacs_AIRS_General/AUG30_2012/cov_10_97_97.mat';
driver.oem.lambda_qst       = [0.1 0.1   0.2 0.2  0.3 0.4 0.5 0.6  1 1];
driver.oem.lambda_Q1        = [ones(1,10) ones(1,30)*1.5 ones(1,57)*1]*100;
driver.oem.lambda_temp      = [ones(1,10)*2 ones(1,30)*1.5 ones(1,57)*1]*50;

driver.oem.lambda_qst       = ones(1,10) * 0.1;
driver.oem.lambda_Q1        = [ones(1,10)*1 ones(1,30)*1.0 ones(1,57)*1]*1;
driver.oem.lambda_temp      = [ones(1,10)*1 ones(1,30)*1.0 ones(1,57)*1]*1;

%% new 9_97_97_97 stuff
driver.jacobian.filename    = '../Aux_jacs_IASI_General/AUG30_2012/all_kcarta_jacsB1_9_97_97_97.mat';
driver.jacobian.filename    = '../Cluster_IASI_General/strow_fake36_Aug2012_IASI_sergio_jac_9_97_97_97fakeB1.mat';
driver.jacobian.filename    = '../Cluster_IASI_General/avg_strow_fake36_Aug2012_IASI_sergio_jac_9_97_97_97fakeB1.mat';
driver.oem.apriori_filename = '../Aux_jacs_AIRS_General/AUG30_2012/apriori300.mat';
driver.oem.cov_filename     = '../Aux_jacs_AIRS_General/AUG30_2012/cov_9_97_97_97.mat';
driver.jacobian.numQlays    = 2;                         %% in addition to water, we are adding on HDO
driver.jacobian.Q2jacindex  = 1:97;
driver.jacobian.qstnames = {'CO2trop' 'CO2strat' 'O3trop' 'O3strat' 'N2O' 'CO' 'CH4' 'CFC11' 'stemp'};
driver.jacobian.qstYesOrNo = [1       1          1        1         1     1     1    1         1];
driver.oem.lambda_qst       = [0.1 0.1   0.1 0.1   0.3 0.4 0.5 0.6  0.1]*1;
driver.oem.lambda_Q1        = [ones(1,10) ones(1,30)*1.5 ones(1,57)*1]*100;
driver.oem.lambda_temp      = [ones(1,10)*2 ones(1,30)*1.5 ones(1,57)*1]*100;

driver.oem.lambda_qst       = [0.1 0.1   0.2 0.2  0.3 0.4 0.5 0.6  1]*0.10*10;
driver.oem.lambda_Q1        = [ones(1,10) ones(1,30)*1.5 ones(1,57)*1]*0.10*10;
driver.oem.lambda_temp      = [ones(1,10)*2 ones(1,30)*1.5 ones(1,57)*1]*10;

driver.oem.diag_only        = 0;
driver.oem.lambda           = 0.01;
driver.oem.lambda_qst       = 0.1;
driver.oem.lambda_Q1     = 0.1;
driver.oem.lambda_temp      = 0.1;

driver.oem.lambda_Q1     = 10;
driver.oem.lambda_temp      = 1;
driver.oem.lambda_Q1     = 50;
driver.oem.lambda_temp      = 10;
driver.oem.lambda_Q1     = 100;
driver.oem.lambda_temp      = 50;

driver.oem.lambda_Q2        = driver.oem.lambda_Q1/10;
%}

%%%%%%%%%%%%%%%%%%%%%%%%%
driver.iibin = ix;
driver.filename = ['../Output/testx_' int2str(driver.iibin)];

