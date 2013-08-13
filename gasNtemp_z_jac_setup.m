function [driver,jac] = gasNtemp_z_jac_setup(driver0)

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

Q1jacindex = zeros(1,numlays);
Q1jacindex(driver.jacobian.Q1jacindex) = 1;
Q1jacindex = logical(Q1jacindex);

if driver.jacobian.numQlays > 1
  for ii = 2 : driver.jacobian.numQlays
    junk = ['Q' num2str(ii) 'jacindex = zeros(1,numlays);']; eval(junk);
    junk = ['Q' num2str(ii) 'jacindex(driver.jacobian.Q' num2str(ii) 'jacindex) = 1;']; eval(junk);
    junk = ['Q' num2str(ii) 'jacindex = logical(Q' num2str(ii) 'jacindex);']; eval(junk);
  end
end

tjacindex = zeros(1,numlays);
tjacindex(driver.jacobian.tjacindex) = 1;
tjacindex = logical(tjacindex);

%% driver.jacindex = [qstjacindex Q1jacindex tjacindex];
driver.jacindex = [];
for ii = 1 : driver.jacobian.numQlays
  junk = ['driver.jacindex = [driver.jacindex Q' num2str(ii) 'jacindex];']; eval(junk);
end
driver.jacindex = logical([qstjacindex driver.jacindex tjacindex]);

%% tracegases/ST + T(z) + [Q1(z) + Q2(Z) + .. QN(z)];
allind = length(qstjacindex) + numlays + numlays*numQprofiles;   
allind = 1 : allind;

goodind = driver.jacindex.*allind;
iqst    = find( goodind > 0 & goodind <= length(qstjacindex));
itemp   = find( goodind > 0 & goodind >= length(qstjacindex) + (driver.jacobian.numQlays*numlays+1));
for ii = 1 : driver.jacobian.numQlays
  aa = length(qstjacindex) + 1;     %% this takes care of qst
  bb = aa + (ii-0)*numlays-1;       %% want to go to end,   hence ii-0
  aa = aa + (ii-1)*numlays;         %% want to go to start, hence ii-1
  junk = ['iQ' num2str(ii) ' = find( goodind > 0 & goodind >= ' num2str(aa) ' & goodind <= ' num2str(bb) ');']; eval(junk)
end

% Find iqst, iwater, itemp for new compressed matrices (done in retrieval.m)
inew = [];
for ii = 1 : driver.jacobian.numQlays
  junk = ['inew = [inew iQ' num2str(ii) '];']; eval(junk);
end
inew = [iqst inew itemp];
[c niqst ib]   = intersect(inew,iqst);
[c nitemp ib]  = intersect(inew,itemp);
for ii = 1 : driver.jacobian.numQlays
  junk = ['[c niQ' num2str(ii) ' ib] = intersect(inew,iQ' num2str(ii) ');']; eval(junk);
end

%%%% 
if findstr(driver.jacobian.filename,'IASI')
  jac1  = squeeze(M_TS_jac_all1(driver.iibin,:,:));
  jac2  = squeeze(M_TS_jac_all2(driver.iibin,:,:));
  jac = [jac1; jac2];
else 
  jac  = squeeze(M_TS_jac_all(driver.iibin,:,:));
end

driver.jacobian.iqst   = niqst;
for ii = 1 : driver.jacobian.numQlays
  junk = ['driver.jacobian.iQ' num2str(ii) ' = niQ' num2str(ii) ';']; eval(junk);
end
driver.jacobian.itemp  = nitemp;

driver.jacobian.gasid   = goodind(iqst);
for ii = 1 : driver.jacobian.numQlays
  junk = ['driver.jacobian.iQ' num2str(ii) 'id = goodind(iQ' num2str(ii) ');']; eval(junk);
end
driver.jacobian.tempid  = goodind(itemp);

%%%%
% Need to avoid these normalizations if just doing LLS and no profile fits
iDoRenorm = -1;   %% skip the renorms
%iDoRenorm = +1;   %% do the renorms
qrenorm0 = qrenorm;
if ((length(iQ1) > 1 | length(itemp) > 1)) & (iDoRenorm > 0)
  disp('RENORM')
  do_the_renorm
end
driver.qrenorm  = qrenorm;
