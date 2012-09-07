function [jac] = get_jacs0(driver);

%% this is very similar to get_jac.m, except it is used when we adjust the rateset
%% eg see Cluster_AIRS/run_retrieval_cluster.m (driver.rateset.adjust = true)
   
% Load Jacobians, assumes they are named M_TS_jac_all
% Also loads qrenorm
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

qstjacindex_names = driver.jacobian.qstnames;    %% this sets the names
qstjacindexX      = driver.jacobian.qstYesOrNo;  %% this sets whether or not to do the jacs
qstjacindex       = logical(qstjacindexX);       %% turn 0/1 into F/T

numlays = driver.jacobian.numlays;        %% how many layers in jacobians eg 97
numQprofiles = driver.jacobian.numQlays;  %% how many gas profiles to retrieve (>= 1for WV at least)

if numQprofiles < 1
  error('need at least one gas profile');
end

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

tjacindex = ones(1,97);
tjacindex = logical(tjacindex);

%% driver.jacindex = [qstjacindex(1:6) wvjacindex tjacindex];
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

%%%% 
if findstr(driver.jacobian.filename,'IASI')
  jac1  = squeeze(M_TS_jac_all1(driver.iibin,:,:));
  jac2  = squeeze(M_TS_jac_all2(driver.iibin,:,:));
  jac = [jac1; jac2];
else 
  jac  = squeeze(M_TS_jac_all(driver.iibin,:,:));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Don't need goodind(iwater), iwater will do..., etc.
for ii = 1 : driver.jacobian.numQlays
  junk = ['iQ' num2str(ii) 'sum = sum(jac(:,goodind(iQ' num2str(ii) ')),2);']; eval(junk);
end
tsum = sum(jac(:,goodind(itemp)),2);


for ii = 1 : driver.jacobian.numQlays
  junk = ['iQ' num2str(ii) 'sum_max = max(abs(iQ' num2str(ii) 'sum));']; eval(junk);
end
tsum_max = max(abs(tsum));

for i=1:length(qstjacindex)
   qsum_max(i) = max(abs(jac(:,i)));
end

% Pick temperature as the standard
for ii = 1 : driver.jacobian.numQlays
  junk = ['q' num2str(ii) '_mult = tsum_max/iQ' num2str(ii) 'sum_max;']; eval(junk);
end
for i=1:length(qstjacindex)
   q_mult(i) = tsum_max/qsum_max(i);
end

% Now apply to Jacobian and modify qrenorm 
for i=1:length(qstjacindex)
   jac(:,i)   = jac(:,i).*q_mult(i);
   qrenorm(i) = qrenorm(i).*q_mult(i);
end
for ii = 1 : driver.jacobian.numQlays
  aa = length(qstjacindex) + 1;     %% this takes care of qst
  bb = aa + (ii-0)*numlays-1;       %% want to go to end,   hence ii-0
  aa = aa + (ii-1)*numlays;         %% want to go to start, hence ii-1
  ind = aa:bb;
  junk = ['themult = q' num2str(ii) '_mult;']; eval(junk);
  jac(:,ind) = jac(:,ind).*themult;
  qrenorm(ind) = qrenorm(ind).*themult;
end
