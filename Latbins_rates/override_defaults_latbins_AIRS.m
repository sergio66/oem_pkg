function driver = override_defaults(driver,ix);

% can use the package eg as
% for ix=1:36                                  
%   ix
%   clear driver; run_retrieval; end

load /strowdata1/shared/sergio/MATLABCODE/TMP_RATES_Fit_pkg/Cluster/strow_stmNov2011_dobs20
driver.jacobian.chanset = dobs20.jacobian.chanset;

%% these are from strow_stmNov2011_dobs20.mat, where Larrabee says he gets
%%    good CH4 rates

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

driver.rateset.datafile = 'fake_transcom_cloudrates.mat';
driver.rateset.datafile = 'fitout_8year_v4_robust.mat';
driver.rateset.datafile = 'fitout_8year_v4_robust_freqcal.mat';
driver.rateset.ncfile   = 'all_lagcor.mat';
driver.rateset.ocb_set  = 'obs';
driver.rateset.adjust   = 0;

%%% larrabee used these in the file he gave me, but I get bad trace gas rates
%%% they give decent WV/T rates
%driver.oem.cov_filename     = 'cov_lls.mat'
driver.oem.diag_only        = 0;
driver.oem.lambda           = 1;
driver.oem.lambda_qst       = 0.1;
driver.oem.lambda_Q1        = 10;
driver.oem.lambda_temp      = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% ABOVE IS STROW, BELOW IS MINE %%%%%%%%%%%
%%%% ABOVE IS STROW, BELOW IS MINE %%%%%%%%%%%

driver.oem.lambda_Q1        = 1;
driver.oem.lambda_temp      = 0.5;

%%%% not bad for all (trace gases and T,WV)
driver.oem.diag_only        = 0;
driver.oem.lambda           = 0.01;
driver.oem.lambda_qst       = 0.1;
driver.oem.lambda_Q1        = 10;
driver.oem.lambda_temp      = 1;
%%%% not bad for all (trace gases and T,WV)

%driver.oem.lambda           = 0.01;
%driver.oem.lambda_temp      = 150;
%driver.oem.lambda_Q1        = 5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

driver.iibin = ix;
driver.filename = ['../Output/testx_' int2str(driver.iibin)];

driver.rateset.datafile = '/strowdata1/shared/sergio/MATLABCODE/RATES_TARO/MAT/overocean_gsx_1dayV1_era2_lays_spanday01_avgL1Brates_robust_June12_2012_span_09_2002_09_2004.mat';   %% 09/01/2002 - 08/31/2004 TO COMPARE AGAINST RSS AMSU
driver.rateset.datafile = '/home/sergio/MATLABCODE/RATES_TARO/MAT/overocean_gsx_1day_clr_era_lays_spanday01_day_avgL1Brates_robust_June29_2012_span_09_2002_08_2011.mat';
driver.rateset.datafile = 'fitout_8year_v4_robust.mat';
driver.rateset.datafile = '/strowdata1/shared/sergio/MATLABCODE/RATES_TARO/MAT/overocean_gsx_1dayV1_era2_lays_spanday01_avgL1Brates_robust_Mar26_2012_span_08_2002_08_2010.mat';   %% 09/01/2002 - 08/31/2010 ERA NIGHT GREAT! as it duplicates LLS fitout_8year_v4_robust.mat
driver.rateset.datafile = '/home/sergio/MATLABCODE/RATES_TARO/MAT/overocean_gsx_1day_clr_era_lays_spanday01_avgL1Brates_robust_June29_2012_span_09_2002_08_2011.mat';

driver.rateset.datafile = 'fitout_8year_v4_robust_freqcal.mat';

driver.rateset.datafile = '/home/sergio/MATLABCODE/RATES_TARO/MAT/overocean_gsx_1day_clr_era_lays_spanday01_avgL1Brates_robust_June29_2012_span_12_2003_11_2011.mat';   %% uses robust, throws out outliers
driver.rateset.datafile = '/home/sergio/MATLABCODE/RATES_TARO/MAT/xoverocean_gsx_1day_clr_era_lays_spanday01_avgL1Brates_robust_June29_2012_span_12_2003_11_2011.mat';   %% uses robust, keeps outliers

driver.rateset.datafile = '/home/sergio/MATLABCODE/RATES_TARO/MAT/overocean_gsx_1day_clr_era_lays_spanday01_avgL1Brates_robust_Aug01_2012_span_01_2003_12_2010.mat';   %% uses robust, throws out outliers
driver.rateset.datafile = '/home/sergio/MATLABCODE/RATES_TARO/MAT/xoverocean_gsx_1day_clr_era_lays_spanday01_avgL1Brates_robust_Aug01_2012_span_01_2003_12_2010.mat';   %% uses robust, keeps outliers

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% orig 6_97_97 stuff
driver.jacobian.filename = 'M_TS_jac_all.mat';
driver.jacobian.qstnames = {'CO2' 'O3' 'N2O' 'CH4' 'CFC11' 'stemp'};
driver.jacobian.qstYesOrNo = [1     1   1     1      1       1];
driver.jacobian.numQlays   = 1;   %% have N lays of temp jacs; how many gases we want to retrieve profiles
                                  %% must be at least 1 (for water)
driver.jacobian.numlays     = 97;
driver.oem.lambda_qst       = [0.1 0.2 0.3 0.4 0.5 0.6]*0.10*10;
driver.oem.lambda_Q1        = [ones(1,10) ones(1,30)*1.5 ones(1,57)*1]*0.10*10;
driver.oem.lambda_temp      = [ones(1,10)*2 ones(1,30)*1.5 ones(1,57)*1]*10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% new 10_97_97 stuff
driver.jacobian.filename = '../Aux_jacs_AIRS_General/AUG30_2012/all_kcarta_jacs_10_97_97.mat';
driver.jacobian.qstnames = {'CO2trop' 'CO2strat' 'O3trop' 'O3strat' 'N2O' 'CO' 'CH4' 'CFC11' 'HDO' 'stemp'};
driver.jacobian.qstYesOrNo = [1       1          1        1         1     1     1    1       1     1];
driver.jacobian.numQlays   = 1;   %% have N lays of temp jacs; how many gases we want to retrieve profiles
                                  %% must be at least 1 (for water)
driver.jacobian.numlays    = 97;
driver.oem.apriori_filename = '../Aux_jacs_AIRS_General/AUG30_2012/apriori204.mat';
driver.oem.cov_filename     = '../Aux_jacs_AIRS_General/AUG30_2012/cov_10_97_97.mat';
driver.oem.lambda_qst       = [0.1 0.1   1 1   0.3 0.4 0.5 0.6  0.1 0.1]*10;
driver.oem.lambda_Q1        = [ones(1,10) ones(1,30)*1.5 ones(1,57)*1]*5;
driver.oem.lambda_temp      = [ones(1,10)*2 ones(1,30)*1.5 ones(1,57)*1]*5;
driver.oem.lambda       = 1

%% new 9_97_97_97 stuff
driver.jacobian.filename    = '../Aux_jacs_AIRS_General/AUG30_2012/all_kcarta_jacs_9_97_97_97.mat';
driver.oem.apriori_filename = '../Aux_jacs_AIRS_General/AUG30_2012/apriori300.mat';
driver.oem.cov_filename     = '../Aux_jacs_AIRS_General/AUG30_2012/cov_9_97_97_97.mat';
driver.rateset.datafile = 'fitout_8year_v4_robust_freqcal.mat';
driver.jacobian.numQlays    = 2;                         %% in addition to water, we are adding on HDO
driver.jacobian.Q2jacindex  = 1:97;
driver.jacobian.qstnames = {'CO2trop' 'CO2strat' 'O3trop' 'O3strat' 'N2O' 'CO' 'CH4' 'CFC11' 'stemp'};
driver.jacobian.qstYesOrNo = [1       1          1        1         1     1     1    1         1];
driver.oem.lambda_qst       = [0.1 0.1   1 1   0.3 0.4 0.5 0.6  0.1]*10;
driver.oem.lambda_Q1        = [ones(1,10) ones(1,30)*1.5 ones(1,57)*1]*1;
driver.oem.lambda_Q2        = driver.oem.lambda_Q1;
driver.oem.lambda_temp      = [ones(1,10)*2 ones(1,30)*1.5 ones(1,57)*1]*1;

driver.rateset.datafile = '/home/sergio/MATLABCODE/RATES_TARO/MAT/xoverocean_gsx_1day_clr_era_lays_spanday01_avgL1Brates_robust_Aug01_2012_span_01_2003_12_2010.mat';   %% uses robust, keeps outliers
driver.oem.lambda_qst       = [ones(1,9)*0.1];
driver.oem.lambda_Q1        = [ones(1,10) ones(1,30)*1.5 ones(1,57)*1]*1;
driver.oem.lambda_Q2        = driver.oem.lambda_Q1;
driver.oem.lambda_temp      = [ones(1,10)*2 ones(1,30)*1.5 ones(1,57)*1]*1;
%%%%%%%%%%%%%%%%%%%%%%%%%

driver.rateset.datafile = '/home/sergio/MATLABCODE/RATES_TARO/MAT/overocean_gsx_1day_clr_era_lays_spanday01_avgL1Brates_robust_Aug28_2012_span_07_2007_07_2012.mat';

fprintf(1,'%s \n',driver.rateset.datafile);
xstartup
g = dogoodchan;
driver.jacobian.chanset = g;