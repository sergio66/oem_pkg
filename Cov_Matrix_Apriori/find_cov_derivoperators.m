function dacov = find_cov_derivoperators(driver);

%% uses the smoothing operators

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

zero   = driver.water_smooth_operator_0;
first  = driver.water_smooth_operator_1;
second = driver.water_smooth_operator_2;
kk = length(wvjacindex);
gas1COV = zero*derivative_operators(kk,0) + ...
     first*derivative_operators(kk,1) + second*derivative_operators(kk,2);
% $$$ gas1COV = gas1COV/driver.qrenorm_wv;
%imagesc(gas1COV); colorbar; pause

zero   = driver.temp_smooth_operator_0;
first  = driver.temp_smooth_operator_1;
second = driver.temp_smooth_operator_2;
kk = length(tjacindex);
tempCOV = zero*derivative_operators(kk,0) + ...
     first*derivative_operators(kk,1) + second*derivative_operators(kk,2);
tempCOV = first*derivative_operators(kk,1) + second*derivative_operators(kk,2);
% $$$ tempCOV = tempCOV/driver.qrenorm_t;
%imagesc(tempCOV); colorbar; pause

ix = ones(1,sum(qstjacindex)); 
% $$$ ix = ix./driver.qrenorm(1:6);
qCOV = diag(ix);
%imagesc(qCOV); colorbar; pause

lala = length(qCOV) + length(gas1COV) + length(tempCOV);
dacov = zeros(lala,lala);
ix = 1:length(qCOV); dacov(ix,ix) = qCOV;
iy = 1:length(gas1COV); ix = length(qCOV)+iy;                 
  dacov(ix,ix) = gas1COV;
iy = 1:length(tempCOV); ix = length(qCOV)+length(gas1COV)+iy; 
  dacov(ix,ix) = tempCOV;
%imagesc(log10(abs(dacov))); colorbar; pause
disp('in find_cov_derivooperator')
keyboard

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

%%% now do dacov = pinv(dacov!!!!) since Steck want a Regularization Matrix,
%%% while I doo pinv(cov)
dacov = pinv(dacov);    %% *******************************

figure(4); clf; pcolor(log10(abs(pinv(dacov)))); colorbar; 
  title('pinv(cov matrix)')
figure(5); clf; plot(diag(pinv(dacov)));
  title('diag(pinv(cov matrix))')

figure(6); clf; pcolor(log10(abs(dacov))); colorbar; 
  %%title(['damping = ' num2str(covdamp)]); pause(0.1)
figure(7); clf; plot(diag(dacov)); pause(0.1); title('diag(cov)')

