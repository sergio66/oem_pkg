function [iUseRetrParam,iUseChan] = get_sequential_indices_params(airsL2_chans_set,airsL2_list_set,airsL2_IDs_set,fairs0,chanIDairs,fuse,inds,xbSave,driver,iSequential,iNXYZLay,iFitTraceGas);

if nargin == 10
  iNXYZLay = 6;
  iFitTraceGas = -1;
elseif nargin == 11
  iFitTraceGas = -1;
end

if iSequential == -1 
  %% nothing to do, use all chans, retrieve all geophys params
  iUseRetrParam  = 1 : length(xbSave);
  iUseChan       = 1 : length(inds);

elseif iSequential == 150
  %% 15 um T(z) and STEMP
  if  iFitTraceGas == -1
    iUseRetrParam = [6 driver.jacobian.temp_i];
  else
    iUseRetrParam = [1 6 driver.jacobian.temp_i];
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%
  %% this uses fuse ~ 500 chans
  iUseChans15 = find(fuse <= 800);
  iUseWindow  = find((fuse > 820 & fuse <= 960) | (fuse > 1225 & fuse <= 1235));
  iUseChan    = union(iUseChans15,iUseWindow);

  %%%%%%%%%%%%%%%%%%%%%%%%%
  %% this uses fairs0 ~ 2645 chans
  iWVBand = find(fairs0 >= 1300);

  iA = airsL2_IDs_set.tz;
  iB = airsL2_IDs_set.stemp;
  iC = union(union(airsL2_IDs_set.cloudphase,airsL2_IDs_set.cirrus),airsL2_IDs_set.emiss_lw);
  iB = union(iB,iC);
  %iB = union(iB,iUseChans15);

  iA = setdiff(iA,iWVBand);
  iB = setdiff(iB,iWVBand);

  %[~,~,iUseChan] = intersect(inds,union(iA,iB));
  [~,iUseChan,~] = intersect(inds,union(iA,iB));

elseif iSequential == 210 %% = 150 + 60
  %% [STEMP], [all 15 um T(z)], [WV(z) lowest only 4 layers]

  if iFitTraceGas == -1
    iUseRetrParam = [6 driver.jacobian.water_i(length(driver.jacobian.water_i)-3:length(driver.jacobian.water_i)) driver.jacobian.temp_i];
  else
    iUseRetrParam = [1 6 driver.jacobian.water_i(length(driver.jacobian.water_i)-3:length(driver.jacobian.water_i)) driver.jacobian.temp_i];
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%
  %% this uses fuse ~ 500 chans
  iUseChans15 = find(fuse <= 800);
  iUseWindow  = find((fuse > 820 & fuse <= 960) | (fuse > 1225 & fuse <= 1235));
  iUseWindow  = find((fuse > 820 & fuse <= 960));
  iUseChan    = union(iUseChans15,iUseWindow);

  %%%%%%%%%%%%%%%%%%%%%%%%%
  %% this uses fairs0 ~ 2645 chans
  iA = airsL2_IDs_set.tz;
  iB = airsL2_IDs_set.stemp;
  iC = union(union(airsL2_IDs_set.cloudphase,airsL2_IDs_set.cirrus),airsL2_IDs_set.emiss_lw);
  iB = union(iB,iC);
  %[~,~,iUseChan] = intersect(inds,union(iA,iB));
  [~,iUseChan,~] = intersect(inds,union(iA,iB));

elseif iSequential == 214  %% = 150 + 60
  %% [STEMP] + [15 um T(z) and WV(z) lowest iNXYZLay layers]

  iUseRetrParam = [6 driver.jacobian.water_i(length(driver.jacobian.water_i)-(iNXYZLay-1):length(driver.jacobian.water_i)) driver.jacobian.temp_i(length(driver.jacobian.temp_i)-(iNXYZLay-1):length(driver.jacobian.temp_i))];

  %%%%%%%%%%%%%%%%%%%%%%%%%
  %% this uses fuse ~ 500 chans
  iUseChan = find((fuse > 780 & fuse <= 960) | (fuse > 1225 & fuse <= 1235));
  iUseChan = find((fuse > 780 & fuse <= 960));

elseif iSequential == 100
  %% 10 um O3(z)
  iUseRetrParam = [driver.jacobian.ozone_i];

  %%%%%%%%%%%%%%%%%%%%%%%%%
  %% this uses fuse ~ 500 chans
  iUse = find(fuse > 1000 & fuse < 1080);
  iUseChan = iUse;

  %%%%%%%%%%%%%%%%%%%%%%%%%
  %% this uses fairs0 ~ 2645 chans
  iA = airsL2_IDs_set.ozone;
  %[~,~,iUseChan] = intersect(inds,iA);
  [~,iUseChan,~] = intersect(inds,iA);

elseif iSequential == 60
  %% 6 um WV(z)
  iUseRetrParam = [driver.jacobian.water_i];

  %%%%%%%%%%%%%%%%%%%%%%%%%
  %% this uses fuse ~ 500 chans
  iUseWindow = find((fuse > 820 & fuse <= 960) | (fuse > 1090 & fuse <= 1280));
  iUseChans6 = find(fuse > 1380 & fuse <= 1700);
  iUseChanX   = union(iUseChans6,iUseWindow);

  iUseWindow = find((fuse > 820 & fuse <= 960));
  iUseChans6 = find(fuse > 1380 & fuse <= 1700);
  iUseChanX   = union(iUseChans6,iUseWindow);

  iUseChan = iUseChanX;
  %%%%%%%%%%%%%%%%%%%%%%%%%
  %% this uses fairs0 ~ 2645 chans
  iA = airsL2_IDs_set.wvz;
  iB = airsL2_IDs_set.cloudclearing;

  %% code uptil June 2023
  %[~,~,iUseChan] = intersect(inds,union(iA,iB));
  [~,iUseChan,~] = intersect(inds,union(iA,iB));

  %%[~,~,iUseChan] = intersect(inds,iA);
  [~,iUseChan,~] = intersect(inds,iA);
  %iUseChan = union(iUseChan,iUseChanX);
end

plot(fuse(iUseChan),'.-'); title(['Chosen Chans at iSequential = ' num2str(iSequential)])
pause(0.1)
