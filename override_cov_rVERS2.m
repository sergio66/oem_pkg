disp('override_cov_r.m')

figure(2); clf; imagesc(r); colorbar
%ret

figure(9); clf; plot(r(1,:),'bo-'); hold
ALLcov = load('/home/sergio/MATLABCODE/RATES_TARO/MAT/overocean_gsx_1day_clr_era_lays_spanday01_fatsummary_Nov02_2012_ALLcov.mat');

%keyboard
cov_temp = squeeze(ALLcov.ALLcov(driver.iibin,:,:));
r=cov_temp*1.0e+6; 
%r(1:6) = cov_temp(1:6)*1.0e+6;
%r(7:103)=cov_temp(7:103)*1.0e+6; 
r(104:200)=cov_temp(104:200)*5.0e+5; 

%figure(5); clf; imagesc(r./r0000); colorbar
%figure(9); clf; imagesc(r); colorbar
%ret


