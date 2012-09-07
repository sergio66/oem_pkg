function dacov = find_cov_damp(driver);

% grows a exp(-((i-j)/damp)^2) for W,T and diag for Q

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

numQjacs  = sum(qstjacindex);
numWVjacs = length(wvjacindex);
numTjacs  = length(tjacindex);

covdamp = driver.covdamp

keyboard

ix = 1 : length(wvjacindex);
ixA = ones(length(wvjacindex),1) * ix;   ixA = ixA/0.01;
ixB = ixA';
gas1COV = ixA-ixB;
gas1COV = gas1COV/covdamp; 
gas1COV = exp(-gas1COV.*gas1COV);

ix = 1 : length(tjacindex);          
ixA = ones(length(tjacindex),1) * ix;    ixA = ixA/0.1;
ixB = ixA';
tempCOV = ixA-ixB;
tempCOV = tempCOV/covdamp; 
tempCOV = exp(-tempCOV.*tempCOV);

ix = ones(1,sum(qstjacindex));     ix = ix./driver.qrenorm(1:6);
qCOV = diag(ix);

lala = length(qCOV) + length(gas1COV) + length(tempCOV);
dacov = zeros(lala,lala);
ix = 1:length(qCOV); dacov(ix,ix) = qCOV;
iy = 1:length(gas1COV); ix = length(qCOV)+iy;                 
  dacov(ix,ix) = gas1COV;
iy = 1:length(tempCOV); ix = length(qCOV)+length(gas1COV)+iy; 
  dacov(ix,ix) = tempCOV;


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

figure(6); clf; pcolor(dacov); colorbar; 
  title(['damping = ' num2str(covdamp)]); pause(0.1)
figure(7); clf; plot(diag(dacov)); pause(0.1); title('diag(cov)')

