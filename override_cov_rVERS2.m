disp('override_cov_r.m')

r0000 = r;
figure(2); clf; imagesc(r); colorbar
%ret

factor = 10
factor = 1
rjunk = zeros(1,200);
figure(9); clf; plot(r(1,:),'bo-'); hold
ALLcov = load('/home/sergio/MATLABCODE/RATES_TARO/MAT/overocean_gsx_1day_clr_era_lays_spanday01_fatsummary_Nov02_2012_ALLcov.mat');

%keyboard
boo = squeeze(ALLcov.ALLcov(driver.iibin,:,:));
r = boo*10e10;

figure(5); clf; imagesc(r./r0000); colorbar
figure(9); clf; imagesc(r); colorbar
ret


