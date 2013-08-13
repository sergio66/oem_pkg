function inv_sap0 = geophysical_covariance(driver);

%% sets up the 200x200 geophysical covariance matrix

l_c = driver.oem.sigma.l_c;

%% humidity = GAS 1
sigma_hum = zeros(1,driver.jacobian.numlays);
%% strat
  L1 = driver.oem.sigma.hum_strat_TOPLAY;
  L2 = driver.oem.sigma.hum_strat_BOTLAY;
  sigma_hum(L1:L2) = driver.oem.sigma.hum_strat_VALUE;
%% trop
  L1 = driver.oem.sigma.hum_trop_TOPLAY;
  L2 = driver.oem.sigma.hum_trop_BOTLAY;
  sigma_hum(L1:L2) = driver.oem.sigma.hum_trop_VALUE;
daInd = 1 : driver.jacobian.numlays;              %% eg 1:97
daInd = daInd + length(driver.jacobian.qstnames); %% eg (1:97)+6
sigma_hum = sigma_hum ./ driver.qrenorm(daInd);

%% temperature
sigma_temp = zeros(1,driver.jacobian.numlays);
%% strat
  L1 = driver.oem.sigma.temp_strat_TOPLAY;
  L2 = driver.oem.sigma.temp_strat_BOTLAY;
  sigma_temp(L1:L2) = driver.oem.sigma.temp_strat_VALUE;
%% upper trop 
  L1 = driver.oem.sigma.temp_upper_trop_TOPLAY; 
  L2 = driver.oem.sigma.temp_upper_trop_BOTLAY; 
  sigma_temp(L1:L2) = driver.oem.sigma.temp_upper_trop_VALUE; 
%% lower trop is just trop  
  L1 = driver.oem.sigma.temp_trop_TOPLAY;
  L2 = driver.oem.sigma.temp_trop_BOTLAY;
  sigma_temp(L1:L2) = driver.oem.sigma.temp_trop_VALUE;
daInd = 1 : driver.jacobian.numlays;                        %% eg 1:97
offset = driver.jacobian.numlays*driver.jacobian.numQlays;  %% eg 3 gases, each 97 layers
offset = offset + length(driver.jacobian.qstnames) ;        %% eg add on col CO2,O3,N2O,CH4,CFC,stemp
daInd = daInd + offset; 
sigma_temp = sigma_temp ./ driver.qrenorm(daInd);

%% QST
sigma_qst = zeros(1,length(driver.jacobian.qstnames));
sigma_qst = driver.oem.sigma.qst;
daInd = 1 : length(driver.jacobian.qstnames); 
sigma_qst = sigma_qst ./ driver.qrenorm(daInd);

i1 = length(sigma_qst)+1;
i2 = length(driver.jacobian.qstnames) ;            %% eg add on col CO2,O3,N2O,CH4,CFC,stemp
i2 = i2 + driver.jacobian.numlays*driver.jacobian.numQlays;  %% eg 3 gases, each 97 layers
iINT = i2;
i2 = i2 + driver.jacobian.numlays;                           %% eg T(97)

s_ap = zeros(i2,i2);
junk = sigma_qst.*sigma_qst;
junk = diag(junk);
s_ap(1:i1-1,1:i1-1) = junk;

for i = i1 : i2
   for j = i1 : i2
      var = 0;
      if i <= iINT & j <= iINT
          var = sigma_hum(i-(i1-1))*sigma_hum(j-(i1-1));
      end
      if i > iINT & j > iINT
          var = sigma_temp(i-iINT)*sigma_temp(j-iINT);
      end
      s_ap(i,j)=var*exp(-(real(i) - real(j))^2/(2*l_c^2));
   end
end
s_ap0    = s_ap;
inv_sap0 = inv(s_ap0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iDebug = -1;
if iDebug>0 
  sigma_temp(1:49)=4;  % Stratophere
  sigma_temp(50:97)=.5; % Troposphere
  sigma_hum(1:49)=1.5;       % Stratosphere
  sigma_hum(50:97)=1;    % Troposphere
  s_ap(1:200,1:200)=0;
  l_c = 2.4;
  for i=7:200
     for j=7:200
        var=0;
        if i<=103 && j<=103
          var=sigma_hum(i-6)*sigma_hum(j-6);
        end
        if i>103 && j >103
          var=sigma_temp(i-103)*sigma_temp(j-103);
        end
        s_ap(i,j)=var*exp(-(real(i) - real(j))^2/(2*l_c^2));
     end
  end
  sigma_co2=4;
  sigma_n2o=4;
  sigma_o3=1;
  sigma_ch4=1;
  sigma_cfc11=1;
  sigma_st=1 ;

  s_ap(1,1)=sigma_co2*sigma_co2;
  s_ap(2,2)=sigma_o3*sigma_o3; ;
  s_ap(3,3)=sigma_n2o*sigma_n2o;
  s_ap(4,4)=sigma_ch4*sigma_ch4;
  s_ap(5,5)=sigma_cfc11*sigma_cfc11;
  s_ap(6,6)=sigma_st*sigma_st;
  disp('keyboard in geophysical_covariance.m so you can check sap matrix')
end

