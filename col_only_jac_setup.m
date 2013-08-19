function [driver,jac] = col_only_jac_setup(driver0)

driver = driver0;

load(driver.jacobian.filename)
driver.f = f;

if findstr(driver.jacobian.filename,'IASI')
  disp('IASI so need tons of jacobians!')
  M_TS_jac_all1 = M_TS_jac_all;
  f1            = f;
  newname = driver.jacobian.filename;
  donk = findstr(newname,'B1');
  newname(donk:donk+1) = 'B2';
  load(newname);
  M_TS_jac_all2 = M_TS_jac_all;
  f2            = f;
  f = [f1 f2];
  driver.f = f;
end

% Form Jacobian selection vector
qstjacindex_names = driver.jacobian.qstnames;    %% this sets the names
qstjacindexX      = driver.jacobian.qstYesOrNo;  %% this sets whether or not to do the jacs
qstjacindex       = logical(qstjacindexX);       %% turn 0/1 into F/T

numlays = driver.jacobian.numlays;        %% how many layers in jacobians eg 97
numQprofiles = driver.jacobian.numQlays;  %% how many gas profiles to retrieve (>= 1for WV at least)

if numlays > 0
  error('huh? called this assuming numlays <= 0')
elseif numQprofiles > 0
  error('huh? called this assuming numQprofiles <= 0')
end

%% driver.jacindex = [qstjacindex Q1jacindex tjacindex] --> qstjacindex
driver.jacindex = logical([qstjacindex]);

%% tracegases/ST + T(z) + [Q1(z) + Q2(Z) + .. QN(z)] ---> tracegases/ST
allind = length(qstjacindex);
allind = 1 : allind;

goodind = driver.jacindex.*allind;
iqst    = find( goodind > 0 & goodind <= length(qstjacindex));

% Find iqst, iwater, itemp for new compressed matrices (done in retrieval.m)
inew = [];
%% inew = [iqst inew itemp]; --> iqst
inew = iqst;
[c niqst ib]   = intersect(inew,iqst);

%%%% 
if findstr(driver.jacobian.filename,'IASI')
  jac1  = squeeze(M_TS_jac_all1(driver.iibin,:,:));
  jac2  = squeeze(M_TS_jac_all2(driver.iibin,:,:));
  jac = [jac1; jac2];
else 
  jac  = squeeze(M_TS_jac_all(driver.iibin,:,:));
end

driver.jacobian.iqst   = niqst;
driver.jacobian.iQ1    = 0;
driver.jacobian.itemp  = 0;

driver.jacobian.gasid   = goodind(iqst);
driver.jacobian.iQ1     = [];
driver.jacobian.tempid  = [];

%%%%
% Need to avoid these normalizations if just doing LLS and no profile fits
iDoRenorm = -1;   %% skip the renorms
%iDoRenorm = +1;   %% do the renorms
qrenorm0 = qrenorm;
if (iDoRenorm > 0)
  disp('RENORM')
  do_the_renorm
end
driver.qrenorm  = qrenorm;
