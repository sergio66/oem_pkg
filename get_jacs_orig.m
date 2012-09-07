function [driver,jac] = get_jacs(driver);
   
% Load Jacobians, assumes they are named M_TS_jac_all
% Also loads qrenorm
load(driver.jacobian.filename)
driver.f = f;
if findstr(driver.jacobian.filename,'IASI')
  disp('IASI so need tons of jacobians!')
  M_TS_jac_all1 = M_TS_jac_all;
  f1            = f;
  newname = driver.jacobian.filename;
  donk = findstr(newname,'B1');
  newname(donk:donk+1) = 'B2';
  load(newname);
  M_TS_jac_all2 = M_TS_jac_all;
  f2            = f;
  f = [f1 f2];
  driver.f = f;
end

% Form Jacobian selection vector
qstjacindex = [driver.jacobian.co2    ...
               driver.jacobian.o3     ...
               driver.jacobian.n2o    ...
               driver.jacobian.ch4    ...
               driver.jacobian.cfc11  ...
               driver.jacobian.stemp];

wvjacindex = zeros(1,97);
wvjacindex(driver.jacobian.wvjacindex) = 1;
wvjacindex = logical(wvjacindex);

tjacindex = zeros(1,97);
tjacindex(driver.jacobian.tjacindex) = 1;
tjacindex = logical(tjacindex);

driver.jacindex = [qstjacindex(1:6) wvjacindex tjacindex];

allind  = 1:200;
goodind = driver.jacindex.*allind;
itemp   = find( goodind > 0 & goodind >= 104);
iwater  = find( goodind > 0 & goodind >= 7 & goodind <= 103);
igas    = find( goodind > 0 & goodind <= 6);

if findstr(driver.jacobian.filename,'IASI')
  jac1  = squeeze(M_TS_jac_all1(driver.iibin,:,:));
  jac2  = squeeze(M_TS_jac_all2(driver.iibin,:,:));
  jac = [jac1; jac2];
else 
  jac  = squeeze(M_TS_jac_all(driver.iibin,:,:));
end

%whos f jac
%figure(2); plot(f,jac(:,1:6))

% Need to avoid these normalizations if just doing LLS and no profile fits
if (length(iwater) > 1 | length(itemp) > 1) 
   
   % Don't need goodind(iwater), iwater will do..., etc.
   % quick ugly fix when not retrieval water
   if length(iwater) > 10
   wsum = sum(jac(:,goodind(iwater)),2);
   else
      wsum = 1;
   end
   
   tsum = sum(jac(:,goodind(itemp)),2);

   wsum_max = max(abs(wsum));
   tsum_max = max(abs(tsum));
   for i=1:6
      qsum_max(i) = max(abs(jac(:,i)));
   end

   % Pick temperature as the standard
   w_mult = tsum_max/wsum_max;
   for i=1:6
      q_mult(i) = tsum_max/qsum_max(i);
   end

   % Now apply to Jacobian and modify qrenorm 
   jac(:,7:104) = jac(:,7:104).*w_mult;
   qrenorm(7:104) = qrenorm(7:104).*w_mult;
   for i=1:6
      jac(:,i)   = jac(:,i).*q_mult(i);
      qrenorm(i) = qrenorm(i).*q_mult(i);
   end

end

driver.qrenorm  = qrenorm;

% Find igas, iwater, itemp for new compressed matrices (done in retrieval.m)
inew = [igas iwater itemp];
[c nigas ib]   = intersect(inew,igas);
[c niwater ib] = intersect(inew,iwater);
[c nitemp ib]  = intersect(inew,itemp);

driver.jacobian.igas   = nigas;
driver.jacobian.iwater = niwater;
driver.jacobian.itemp  = nitemp;

driver.jacobian.gasid   = goodind(igas);
driver.jacobian.waterid = goodind(iwater);
driver.jacobian.tempid  = goodind(itemp);

