function [driver,jac] = get_jacs(driver);
   
% Load Jacobians, assumes they are named M_TS_jac_all
% Also loads qrenorm

numlays = driver.jacobian.numlays;        %% how many layers in jacobians eg 97
numQprofiles = driver.jacobian.numQlays;  %% how many gas profiles to retrieve (>= 1for WV at least)

if numQprofiles < 1 | numlays < 1
  disp('numQprofiles < 1 | numlays < 0 ==> code is doing COL Q jacobians and/or stemp')
end

if numQprofiles >= 1 & numlays >= 1
  %% definitely doing at least one gas(z) jacobian, as well as T(z) jacobian
  [driver,jac] = gasNtemp_z_jac_setup(driver);
elseif numQprofiles < 1 & numlays < 1
  %% doing col jacs
  [driver,jac] = col_only_jac_setup(driver);
elseif numQprofiles < 1 & numlays > 1
  %% doing col jacs for gases, stemp; and doing a T(z) jacobian
  error('[driver,jac] = gascol_temp_z_setup(driver); NOT YET DONE')
end
