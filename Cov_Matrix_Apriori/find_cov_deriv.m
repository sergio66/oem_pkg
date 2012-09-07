function dacov = find_cov_deriv(driver)

loader = ['load ' driver.largesummaryfile]; eval(loader)

addpath /asl/matlab/science
datime = squeeze(nanmean(a.rtime_avg,2));
  [yy,mm,dd,hh] = tai2utc(datime);
days = (yy-2002)*365 + (mm-1)*30 + dd;
yymm = yy + (mm-1)/12;    %% so 01/2002 becomes 2002.0, 07/2002 becomes 2002.5 and 12/2002 becomes 2002.91

d16days = days(1):16:days(length(days));

ii = driver.iibin;
plays = squeeze(a.plevs_avg(:,ii,:));
stemp = squeeze(a.stemp_avg(1,ii,:));
temp = squeeze(a.ptemp_avg(:,ii,:));
gas1 = squeeze(a.gas_1_avg(:,ii,:));
gas2 = squeeze(a.gas_2_avg(:,ii,:));
gas3 = squeeze(a.gas_3_avg(:,ii,:));
gas4 = squeeze(a.gas_4_avg(:,ii,:));
gas5 = squeeze(a.gas_5_avg(:,ii,:));
gas6 = squeeze(a.gas_6_avg(:,ii,:));

stempx = interp1(days,stemp,d16days,[],'extrap');
playsx = interp1(days,plays',d16days,[],'extrap'); playsx = playsx';
tempx = interp1(days,temp',d16days,[],'extrap'); tempx = tempx';
gas1x = interp1(days,gas1',d16days,[],'extrap'); gas1x = gas1x';
gas2x = interp1(days,gas2',d16days,[],'extrap'); gas2x = gas2x';
gas3x = interp1(days,gas3',d16days,[],'extrap'); gas3x = gas3x';
gas4x = interp1(days,gas4',d16days,[],'extrap'); gas4x = gas4x';
gas5x = interp1(days,gas5',d16days,[],'extrap'); gas5x = gas5x';
gas6x = interp1(days,gas6',d16days,[],'extrap'); gas6x = gas6x';

lala = (nanmean(gas1x'))'*ones(1,length(d16days)); gas1x = gas1x./lala;
lala = (nanmean(gas2x'))'*ones(1,length(d16days)); gas2x = gas2x./lala;
lala = (nanmean(gas3x'))'*ones(1,length(d16days)); gas3x = gas3x./lala;
lala = (nanmean(gas4x'))'*ones(1,length(d16days)); gas4x = gas4x./lala;
lala = (nanmean(gas5x'))'*ones(1,length(d16days)); gas5x = gas5x./lala;
lala = (nanmean(gas6x'))'*ones(1,length(d16days)); gas6x = gas6x./lala;

%% the derivative = blah/16 day^-1 = blah/16/365 yr^-1
stempx = diff(stempx'); stempx = stempx'/365/16;
tempx = diff(tempx');   tempx = tempx'/365/16;
gas1x = diff(gas1x');   gas1x = gas1x'/365/16;
gas2x = diff(gas2x');   gas2x = gas2x'/365/16;
gas3x = diff(gas3x');   gas3x = gas3x'/365/16;
gas4x = diff(gas4x');   gas4x = gas4x'/365/16;
gas5x = diff(gas5x');   gas5x = gas5x'/365/16;
gas6x = diff(gas6x');   gas6x = gas6x'/365/16;

stemp = stempx; temp = tempx;
gas1 = gas1x; gas2 = gas2x; gas3 = gas3x;
gas4 = gas4x; gas5 = gas5x; gas6 = gas6x;

clear stempx tempx gas*x

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
    gas1x = gas1(ix,:)/0.01;
    gas1COV = [gas1COV; gas1x];
  end
end

tempCOV = [];
for ix = 1 : length(tjacindex)
  if tjacindex(ix)
    txx = temp(ix,:)/0.1;
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
stemp = stemp;

qCOV = [];
if qstjacindex(1)
  qCOV = [qCOV; gas2/driver.qrenorm(1)];
end
if qstjacindex(2)
  qCOV = [qCOV; gas3/driver.qrenorm(2)];
end
if qstjacindex(3)
  qCOV = [qCOV; gas4/driver.qrenorm(3)];
end
if qstjacindex(4)
  qCOV = [qCOV; gas6/driver.qrenorm(4)];
end
if qstjacindex(5)
  qCOV = [qCOV; gas5/driver.qrenorm(5)];  %%% use CO instead of CFC11
end
if qstjacindex(6)
  qCOV = [qCOV; stemp/0.1];
end

dacov = nancov([qCOV; gas1COV; tempCOV]');

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
figure(7); clf; plot(diag(dacov)); pause(0.1); title('diag(cov)')

