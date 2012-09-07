function driver = get_rates_structure();

%% this is the HUGE file that has strcture with profiles, radiances, 
%% times, angles etc
%% by loading this in, the code can generate covariance matrices on the fly
xfile = ['/Users/strow/Desktop/Rates/Data/'];
xfile = [xfile 'overocean_gsx_1dayV1_era2_lays_fatsummary_June14_2011.mat'];
driver.largesummaryfile = xfile;

%% this file contains the rates derived from the ERA profiles
%% as well as the std.devs of the rates. So this can also be used to 
%% generate covariance matrices.
ecmfile = ['/Users/strow/Desktop/Rates/Data/'];
ecmfile = [ecmfile ...
      'overocean_gsx_1dayV1_era2_lays_profilerates_June14_2011_robust.mat'];
driver.ecmfile = ecmfile;

anom_ecmfile = ['/Users/strow/Work/Rates/Geo_fits/'];
anom_ecmfile = [anom_ecmfile ...
 'overocean_gsx_1dayV1_era2_lays_profilerates_June14_2011_robust_anomaly.mat'];
driver.anom_ecmfile = anom_ecmfile;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Which latitude bin
driver.iibin = 4;

%% Choose variables to fit
driver.jacobian.co2   = true;
driver.jacobian.o3    = true;
driver.jacobian.n2o   = true;
driver.jacobian.ch4   = true;
driver.jacobian.cfc11 = true;
driver.jacobian.stemp = true;
%% wv and t have maximum 97 layers (1=toa,97 = gnd)
driver.jacobian.wvjacindex  = 1:97;
driver.jacobian.tjacindex   = 1:97;

%% damping term, if needed (should be positive)
driver.covdamp = 1.0e+4;

%% overall scale factor, if needed (should be positive)
driver.scale = 1e+2;

% which type of covariance matrix to make
% (0) from ECMWF levels --> layers, use cov(anomaly(atmospheric params))
% (1) from ECMWF levels --> layers, use cov(atmospheric params)
% (2) from ECMWF levels --> layers, use cov(d/dt atmospheric params)
% (3) from exp(-((i-j)/d).^2
% (4) from exp(-((i-j)/d).^2 with ECMWF inputs
% (5) from Stecks smoothing operators
driver.cov_choice = 5;

% should we block diagonalize?
driver.block_diagnol = true;

% smoothing operators for water : wgts for 'zero','first', or 'second' derivs?
driver.water_smooth_operator_0 = 0;
driver.water_smooth_operator_1 = 1;
driver.water_smooth_operator_2 = 0;

% smoothing operators for temp : wgts for 'zero','first', or 'second' derivs?
driver.temp_smooth_operator_0 = 0;
driver.temp_smooth_operator_1 = 1;
driver.temp_smooth_operator_2 = 0;

%% this is the renormalization used in the jacs
driver.qrenorm = [[2.2 0.01 1.0 5 1] [0.1] [0.01*ones(1,97)] [0.1*ones(1,97)]];
driver.qrenorm_wv = 0.01;
driver.qrenorm_t  = 0.1;

%driver.qrenorm = ones(size(driver.qrenorm));
%driver.qrenorm_wv = 1.0;
%driver.qrenorm_t  = 1.0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output filename for apriori files
driver.apriorifilename = ['JUNK/apriori_' int2str(driver.iibin)];

%% Output filename for cov files
driver.covfilename = ['JUNK/xcov_latbin_' int2str(driver.iibin) '_covchoice_' int2str(driver.cov_choice)];

