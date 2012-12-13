disp('override_cov_r.m')

figure(2); clf; imagesc(r); colorbar

factor = 10
factor = 1
rjunk = zeros(1,200);
figure(9); clf; plot(r(1,:),'bo-'); hold
co2cov = load('/home/sergio/MATLABCODE/RATES_TARO/MAT/overocean_gsx_1day_clr_era_lays_spanday01_fatsummary_Nov02_2012_co2cov.mat');
co2cov
rjunk(7:47) = co2cov.co2cov(driver.iibin,4);      % WV S
rjunk(48:103) = co2cov.co2cov(driver.iibin,5);    % WV T
rjunk(104:144) = co2cov.co2cov(driver.iibin,1);   % temp S
rjunk(145:200) = co2cov.co2cov(driver.iibin,2);   % temp T
rjunk(6) = co2cov.co2cov(driver.iibin,3);         % stemp
figure(9); clf; plot(rjunk)          

%r(1,:) = r(1,:) + rjunk*factor;
%  plot(r(1,:),'r'); hold off
%r(:,1) = r(:,1) + rjunk'*factor;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

factor = 1
rjunk = zeros(1,200);
figure(9); clf; plot(r(6,:),'bo-'); hold
STcov = load('/home/sergio/MATLABCODE/RATES_TARO/MAT/overocean_gsx_1day_clr_era_lays_spanday01_fatsummary_Nov02_2012_STcov.mat');
STcov
rjunk(7:47) = STcov.STcov(driver.iibin,4);      % WV S
rjunk(48:103) = STcov.STcov(driver.iibin,5);    % WV T
rjunk(104:144) = STcov.STcov(driver.iibin,1);   % temp S
rjunk(145:200) = STcov.STcov(driver.iibin,2);   % temp T
rjunk(6) = STcov.STcov(driver.iibin,3);         % stemp
figure(9); clf; plot(rjunk)          

r(6,:) = r(6,:) + rjunk*factor;
  plot(r(6,:),'r'); hold off
r(:,6) = r(:,6) + rjunk'*factor;

imagesc(r);  colorbar
