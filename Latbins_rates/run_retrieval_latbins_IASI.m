% do_retrieval_hand.m
%
% Script to do OEM and LLS retrieval of radiance rates
% Retrieval parameters set in get_struct
%
% V0.1 Sergio, Larrabee: Sept 11, 2012
%
% can easily do all latbins via
%     for ix = 1 : 36
%       run_retrieval_latbins_IASI
%     end

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
%driver = override_defaults_iasi(driver,ix);  %% this is for LLS oem_pkg
driver = override_defaults_latbins_IASI(driver,ix);

%% cd ../WORKS_Sept1_2012

% Get rate data and Jacobians
driver            = get_rates(driver);
[driver,m_ts_jac] = get_jacs(driver);          %% strow's new renorm
%[driver,m_ts_jac] = get_jacs_NOrenorm(driver); %% no renorm

%  Adjust the rates?
if driver.rateset.adjust
  % m_ts_jac0 = get_jacs_6_97_97(driver);
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

