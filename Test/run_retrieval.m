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

d.rateset.ocb_set      = 'obs';
d.iibin                = 18;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% change some defaults; make sure you correctly choose AIRS vs IASI
%% edit the settings in these files as wanted

d = override_defaults_airs(d,d.iibin);
%d = override_defaults_iasi(d,d.iibin);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get rates and and jacs
d        = get_rates(d);
[d,jacs] = get_jacs(d);

% Adjust the rates?
if d.rateset.adjust
  jacs0 = get_jacs0(d);
  d = adjust_rates(d,jacs0);
end

disp('read in stuff ... doing retrieval ...')

% Do the retrieval 
d = retrieval(d,jacs);

% Save retrieval output
save(d.filename,'-struct','d');

% Close debug file
if d.debug
  writelog('close')
end

