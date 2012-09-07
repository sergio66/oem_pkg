function [dacov] = find_cov_simple_anom(driver);

loader = ['load ' driver.anom_ecmfile]; eval(loader)

%% keyboard

addpath /asl/matlab/science
datime = anomaly.time;   %% time since Jan 1,2002
yymm = datime/365;

ii = driver.iibin;
%plays = squeeze(anomaly.plevs(:,ii,:));
stemp = (anomaly.stemp(ii,:));
temp = squeeze(anomaly.ptemp(:,ii,:));
gas1 = squeeze(anomaly.water(:,ii,:));
gas2 = squeeze(anomaly.co2(:,ii,:));
gas3 = squeeze(anomaly.ozone(:,ii,:));
gas4 = squeeze(anomaly.n2o(:,ii,:));
gas5 = squeeze(anomaly.co(:,ii,:));
gas6 = squeeze(anomaly.ch4(:,ii,:));

%% Determine which Jacobians we need
qstjacindex = [driver.jacobian.co2    ...
               driver.jacobian.o3     ...
               driver.jacobian.n2o    ...
               driver.jacobian.ch4    ...
               driver.jacobian.cfc11  ...
               driver.jacobian.stemp];

wvjacindex = zeros(1,97);
wvjacindex(driver.jacobian.wvjacindex) = 1;

tjacindex = zeros(1,97);
tjacindex(driver.jacobian.tjacindex) = 1;

gas1COV = [];
for ix = 1 : length(wvjacindex)
  if wvjacindex(ix)
    gas1x = gas1(ix,:)/driver.qrenorm_wv;
    gas1COV = [gas1COV; gas1x];
  end
end

tempCOV = [];
for ix = 1 : length(tjacindex)
  if tjacindex(ix)
    txx = temp(ix,:)/driver.qrenorm_t;
    tempCOV = [tempCOV; txx];
  end
end

%% already normalized above so just need to take avg
%gas2 = nansum(gas2(1:95,:))/95; 
%gas3 = nansum(gas3(1:95,:))/95; 
%gas4 = nansum(gas4(1:95,:))/95; 
%gas5 = nansum(gas5(1:95,:))/95; 
%gas6 = nansum(gas6(1:95,:))/95; 
gas2 = nanmean(gas2(1:95,:)); 
gas3 = nanmean(gas3(1:95,:)); 
gas4 = nanmean(gas4(1:95,:)); 
gas5 = nanmean(gas5(1:95,:)); 
gas6 = nanmean(gas6(1:95,:)); 
stemp = stemp';

qCOV = [];
if qstjacindex(1)
  qCOV = [qCOV gas2'/driver.qrenorm(1)];
end
if qstjacindex(2)
  qCOV = [qCOV gas3'/driver.qrenorm(2)];
end
if qstjacindex(3)
  qCOV = [qCOV gas4'/driver.qrenorm(3)];
end
if qstjacindex(4)
  qCOV = [qCOV gas6'/driver.qrenorm(4)];
end
if qstjacindex(5)
  qCOV = [qCOV gas5'/driver.qrenorm(5)];  %%% use CO instead of CFC11
end
if qstjacindex(6)
  qCOV = [qCOV stemp/0.1];
end

dacov = nancov([qCOV'; gas1COV; tempCOV]');

if driver.block_diagnol
  dacov2 = zeros(size(dacov));
  n1 = sum(qstjacindex);
  n2 = sum(wvjacindex);
  n3 = sum(tjacindex);
  [aajunk,bbjunk] = size(dacov);
  if aajunk ~= (n1 + n2 + n3)
    [aajunk n1 n2 n3]
    error('ooops inconsistent!');
  end
  %% block diagonalize the trace gases/ST
  aajunk = 1:n1;
  dacov2(aajunk,aajunk) = dacov(aajunk,aajunk);

  %% block diagonalize the WV
  aajunk = n1+1:n1+n2;
  dacov2(aajunk,aajunk) = dacov(aajunk,aajunk);

  %% block diagonalize the trace gases/ST
  aajunk = n1+n2+1:n1+n2+n3;
  dacov2(aajunk,aajunk) = dacov(aajunk,aajunk);

  dacov = dacov2;
end

figure(6); clf; pcolor(log10(abs(double(dacov)))); colorbar; 
title('covariance matr'); pause(0.1)
figure(7); clf; plot(diag(dacov),'r'); pause(0.1); title('diag(cov)')

