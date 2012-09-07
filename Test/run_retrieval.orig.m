% do_retrieval_hand.m
%
% Script to do OEM and LLS retrieval of radiance rates
% Retrieval parameters set in get_struct
%
% V0.1 Sergio, Larrabee: Aug 1, 2011

%%%% addpath ../                     %% this is my version
%%%% addpath ../Strow_Dec17_2011     %% this is what Strow has
addpath ../                     %% this is my version

% Define default driver structure
d = set_struct;

% Open debug file if desired
if d.debug
  writelog('open');
end;

d.oem.cov_filename     = 'cov_lls.mat';
d.oem.apriori_filename = 'apriori_zero';

d.rateset.ocb_set      = 'obs';
d.iibin                = 20;

% Regularization
d.oem.lambda = 0.1;
switch d.rateset.ocb_set
   case {'obs','cal'}
      d.oem.lambda_gas   = 0.1;
      d.oem.lambda_water = 0.9; 
      d.oem.lambda_temp  = 0.9;
  case 'bias'
      d.oem.lambda_gas   = 0.1;
      d.oem.lambda_water = 5; 
      d.oem.lambda_temp  = 5;
end
      
ig = d.jacobian.chanset;
ig = ig( ig < 2150);
%ig = ig( ig < 2116);
d.jacobian.chanset = ig;

% Remove CO2 from obs, bias
switch d.rateset.ocb_set
case{'obs','bias'}
   d.rateset.adjust        = true;
   d.rateset.adjust_index  = [1];
   d.rateset.adjust_values = [1.986];
   d.jacobian.co2          = false;
case 'cal'
   d.rateset.adjust        = false;
   d.jacobian.co2          = false;
   d.jacobian.o3           = true;
   d.jacobian.n2o          = false;
   d.jacobian.ch4          = false;
   d.jacobian.cfc11        = false;
   d.jacobian.stemp        = true;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% this is for debugging, get rid of it
%% d = override_defaultsX(d);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get rates and and jacs
d        = get_rates(d);
[d,jacs] = get_jacs(d);

% Adjust the rates?
if d.rateset.adjust
  d = adjust_rates(d,jacs);
end

which retrieval
% Do the retrieval 
d = retrieval(d,jacs);

% Save retrieval output
save(d.filename,'-struct','d');

% Close debug file
if d.debug
  writelog('close')
end

