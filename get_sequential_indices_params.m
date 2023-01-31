function [iUseRetrParam,iUseChan] = get_sequential_indices_params(airsL2_chans_set,airsL2_list_set,airsL2_IDs_set,fairs0,chanIDairs,fuse,inds,xbSave,driver,iSequential);

if iSequential == -1 
  %% nothing to do, use all chans, retrieve all geophys params
  iUseRetrParam  = 1 : length(xbSave);
  iUseChan       = 1 : length(inds);

elseif iSequential == 150
  %% 15 um T(z) and STEMP
  iUseRetrParam = [6 driver.jacobian.temp_i];

  iUseChans15 = find(fuse <= 800);
  iUseWindow  = find((fuse > 820 & fuse <= 960) | (fuse > 1225 & fuse <= 1235));
  iUseChan    = union(iUseChans15,iUseWindow);

  iA = airsL2_IDs_set.tz;
  iB = airsL2_IDs_set.stemp;
  iC = union(union(airsL2_IDs_set.cloudphase,airsL2_IDs_set.cirrus),airsL2_IDs_set.emiss_lw);
  iB = union(iB,iC);
  %[~,~,iUseChan] = intersect(inds,union(iA,iB));
  [~,iUseChan,~] = intersect(inds,union(iA,iB));

elseif iSequential == 210 %% = 150 + 60
  %% 15 um T(z) and STEMP and WV(z) lowest 4 layers
  iUseRetrParam = [6 driver.jacobian.water_i(length(driver.jacobian.water_i)-3:length(driver.jacobian.water_i)) driver.jacobian.temp_i];

  iUseChans15 = find(fuse <= 800);
  iUseWindow  = find((fuse > 820 & fuse <= 960) | (fuse > 1225 & fuse <= 1235));
  iUseChan    = union(iUseChans15,iUseWindow);

  iA = airsL2_IDs_set.tz;
  iB = airsL2_IDs_set.stemp;
  iC = union(union(airsL2_IDs_set.cloudphase,airsL2_IDs_set.cirrus),airsL2_IDs_set.emiss_lw);
  iB = union(iB,iC);
  %[~,~,iUseChan] = intersect(inds,union(iA,iB));
  [~,iUseChan,~] = intersect(inds,union(iA,iB));

elseif iSequential == 214  %% = 150 + 60
  %% STEMP, 15 um T(z) and WV(z) lowest 4 layers
  iUseRetrParam = [6 driver.jacobian.water_i(length(driver.jacobian.water_i)-3:length(driver.jacobian.water_i)) driver.jacobian.temp_i(length(driver.jacobian.temp_i)-3:length(driver.jacobian.temp_i))];

  iUseChan = find((fuse > 780 & fuse <= 960) | (fuse > 1225 & fuse <= 1235));

elseif iSequential == 100
  %% 10 um O3(z)
  iUseRetrParam = [driver.jacobian.ozone_i];

  iUse = find(fuse > 1000 & fuse < 1080);
  iUseChan = iUse;

  iA = airsL2_IDs_set.ozone;
  %[~,~,iUseChan] = intersect(inds,iA);
  [~,iUseChan,~] = intersect(inds,iA);

elseif iSequential == 60
  %% 6 um WV(z)
  iUseRetrParam = [driver.jacobian.water_i];

  iUseWindow = find((fuse > 820 & fuse <= 960) | (fuse > 1090 & fuse <= 1280));
  iUseChans6 = find(fuse > 1380 & fuse <= 1700);
  iUseChan   = union(iUseChans6,iUseWindow);

  iA = airsL2_IDs_set.wvz;
  iB = airsL2_IDs_set.cloudclearing;
  %[~,~,iUseChan] = intersect(inds,union(iA,iB));
  [~,iUseChan,~] = intersect(inds,union(iA,iB));

end
