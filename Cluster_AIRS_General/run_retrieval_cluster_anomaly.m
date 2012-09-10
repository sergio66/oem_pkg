% do_retrieval_hand.m
%
% Script to do OEM and LLS retrieval of radiance rates
% Retrieval parameters set in get_struct
%
% V0.1 Sergio, Larrabee: Aug 1, 2011

pwd

addpath ../
addpath ../Test
addpath ../Test/Debug

%% Define default driver structure
driver = set_struct;

%% Open debug file if desired
if driver.debug
  writelog('open');
end;

%% Change some defaults
driver = override_defaults(driver,ix);

% Get rate data and Jacobians
driver            = get_rates(driver);

%%%%% NEW NEW NEW NEW NEW NEW %%% replace rateset with smooth_anomaly
% rx = squeeze(smoothAnomObs(:,iSmoothStep,ix));  %% 380
% rx = squeeze(obs_anomaly(:,iSmoothStep,ix));      %% 3800
rx = smoothed_obs_anomaly(:,iSmoothStep);         %% 3800, smoothed
driver.rateset.rates = rx;
driver.rateset.datafile = smooth_anomaly_file;
%%%%% NEW NEW NEW NEW NEW NEW %%% replace rateset with smooth_anomaly

[driver,m_ts_jac] = get_jacs(driver);          %% strow's new renorm
%[driver,m_ts_jac] = get_jacs_NOrenorm(driver); %% no renorm

%  Adjust the rates?
if driver.rateset.adjust
  m_ts_jac0 = get_jacs0(driver);
  driver = adjust_rates(driver,m_ts_jac0);
end

% Do the retrieval
driver = retrieval(driver,m_ts_jac);

%% Save retrieval output
save(driver.filename,'-struct','driver');

%% Close debug file
if driver.debug
  writelog('close')
end

