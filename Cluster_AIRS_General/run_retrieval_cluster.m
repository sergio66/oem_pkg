% do_retrieval_hand.m
%
% Script to do OEM and LLS retrieval of radiance rates
% Retrieval parameters set in get_struct
%
% V0.1 Sergio, Larrabee: Aug 1, 2011

% pwd

addpath /home/sergio/MATLABCODE/oem_pkg/
addpath /home/sergio/MATLABCODE/oem_pkg/Test/
addpath /home/sergio/MATLABCODE/oem_pkg/Test/Debug

%% Define default driver structure
driver = set_struct;

%% Open debug file if desired
if driver.debug
  writelog('open');
end;

%% Change some defaults
%driver = override_defaults_airs(driver,ix);  %% this is for LLS oem_pkg
driver = override_defaults(driver,ix);

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

