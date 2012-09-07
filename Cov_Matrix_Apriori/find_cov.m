function [dacov,decouple] = find_cov(QSTjacindex,WVjacindex,Tjacindex,latBin)

%keyboard

get_struct

addpath /asl/matlab/science
datime = squeeze(nanmean(a.rtime_avg,2));
  [yy,mm,dd,hh] = tai2utc(datime);
days = (yy-2002)*365 + (mm-1)*30 + dd;
yymm = yy + (mm-1)/12;    %% so 01/2002 becomes 2002.0, 07/2002 becomes 2002.5 and 12/2002 becomes 2002.91

ii = latBin;
plays = squeeze(a.plevs_avg(:,ii,:));
stemp = squeeze(a.stemp_avg(1,ii,:));
temp = squeeze(a.ptemp_avg(:,ii,:));
gas1 = squeeze(a.gas_1_avg(:,ii,:));
gas2 = squeeze(a.gas_2_avg(:,ii,:));
gas3 = squeeze(a.gas_3_avg(:,ii,:));
gas4 = squeeze(a.gas_4_avg(:,ii,:));
gas5 = squeeze(a.gas_5_avg(:,ii,:));
gas6 = squeeze(a.gas_6_avg(:,ii,:));

lala = (nanmean(gas1'))'*ones(1,length(yymm)); gas1 = gas1./lala;
lala = (nanmean(gas2'))'*ones(1,length(yymm)); gas2 = gas2./lala;
lala = (nanmean(gas3'))'*ones(1,length(yymm)); gas3 = gas3./lala;
lala = (nanmean(gas4'))'*ones(1,length(yymm)); gas4 = gas4./lala;
lala = (nanmean(gas5'))'*ones(1,length(yymm)); gas5 = gas5./lala;
lala = (nanmean(gas6'))'*ones(1,length(yymm)); gas6 = gas6./lala;

numQjacs  = sum(QSTjacindex);
numWVjacs = length(WVjacindex);
numTjacs  = length(Tjacindex);

if numWVjacs >= 97
  gas1COV = [];
  for ix = 1 : length(WVjacindex)
    gas1x = gas1(WVjacindex(ix),:);
    gas1COV = [gas1COV; gas1x];
    end
else
  gas1COV = [];
  for ix = 1 : length(WVjacindex)
    boop = WVjacindex(ix);
    gas1x = gas1(WVjacindex(ix),:);
    gas1COV = [gas1COV; gas1x];
    end
  end

if numTjacs >= 97
  tempCOV = [];
  for ix = 1 : length(Tjacindex)
    txx = temp(Tjacindex(ix),:);
    tempCOV = [tempCOV; txx];
    end
else
  tempCOV = [];
  for ix = 1 : length(Tjacindex)
    boop = Tjacindex(ix);
    txx = temp(Tjacindex(ix),:);
    tempCOV = [tempCOV; txx];
    end
  end

%% already normalized above so no need to normalize
gas2 = nansum(gas2(1:95,:)); 
gas3 = nansum(gas3(1:95,:)); 
gas4 = nansum(gas4(1:95,:)); 
gas5 = nansum(gas5(1:95,:)); 
gas6 = nansum(gas6(1:95,:)); 
stemp = stemp';

if numQjacs == 6
  qCOV = [gas2; gas3; gas4; gas6; gas5; stemp]; 
else
  qCOV = [];
  for ix = 1 : length(QSTjacindex)
    if QSTjacindex(ix) == 1
      qCOV = [qCOV; gas2];
      end
    if QSTjacindex(ix) == 1
      qCOV = [qCOV; gas3];
      end
    if QSTjacindex(ix) == 1
      qCOV = [qCOV; gas4];
      end
    if QSTjacindex(ix) == 1
      qCOV = [qCOV; gas6];
      end
    if QSTjacindex(ix) == 1
      qCOV = [qCOV; gas5];
      end
    if QSTjacindex(ix) == 1
      qCOV = [qCOV; stemp];
      end
    end
  end

dacov = nancov([qCOV; gas1COV; tempCOV]');

disp('  cov matrix done from nancov([qCOV; gas1COV; tempCOV])')
disp('  which means q,WV and T are all coupled ie this is NOT block diagnol')
iDecouple = input('  Tweak the cov matrix?? ');
if iDecouple > 0
  iDecouple = input(' --> decouple gases COV from WV/T and from each other? ');
  if iDecouple > 0
    decouple.separategases = +1;
    renormQ = [2.2 0.01 1.0 5 1 0.1];
    ecmwf_rates_gas     = [2.2 0.01  1.0 5 -1.0 0.01]; 
    ecmwf_rates_unc_gas = [0.2 0.001 0.2 2 0.2  0.001]; 
    ecmwf_rates_unc = ecmwf_rates_unc_gas(QSTjacindex)./renormQ(QSTjacindex);
    qCOV = diag(ecmwf_rates_unc.^2);
    dacov(1:length(QSTjacindex),1:length(QSTjacindex)) = qCOV;
    end
  iDecouple = input(' --> decouple WV from T? ');
  if iDecouple > 0
    decouple.separateWVfromT = +1;
    dacov0 = dacov;
    mmQ = length(QSTjacindex);
    mmWV = length(WVjacindex);
    a1 = mmQ+mmWV; a2 = length(dacov);
    dacov(mmQ+1:a1,a1+1:a2) = 0.0;
    dacov(a1+1:a2,mmQ+1:a1) = 0.0;
    end
  iDecouple = input(' --> renormalize so Q,WV,T all same order of mag? ');
  if iDecouple > 0
    decouple.renormQ_WV_T = +1;
    dacov0 = dacov;
    mmQ = length(QSTjacindex);
    mmWV = length(WVjacindex);
    mmT  = length(dacov);
    aaQ = dacov(1:mmQ,1:mmQ);                    aaQ = max(aaQ(:));
    aaWV = dacov(mmQ+1:mmQ+mmWV,mmQ+1:mmQ+mmWV); aaWV = max(aaWV(:));
    aaT = dacov(mmQ+mmWV+1:mmT,mmQ+mmWV+1:mmT);  aaT  = max(aaT(:));
    blah = max([aaQ aaWV aaT]);
    ind = 1: mmQ;         dacov(ind,ind) = dacov(ind,ind) * blah/aaQ/blah;
    ind = mmQ+1:mmQ+mmWV; dacov(ind,ind) = dacov(ind,ind) * blah/aaWV/blah;
    ind = mmQ+mmWV+1:mmT; dacov(ind,ind) = dacov(ind,ind) * blah/aaT/blah;
    end
  iDecouple = input(' --> put blanket off diagnol damping??? ');
  if iDecouple > 0
    dacov0 = dacov;
    covDamp = input('enter damping in exp(-((i-j)/d)^2) : ');
    decouple.covDamp = covDamp;
    ix = 1: length(dacov);
    ixA = ones(length(dacov),1) * ix;
    ixB = ixA';
    tempCOV = ixA-ixB;
    tempCOV = tempCOV/covDamp; 
    tempCOV = exp(-tempCOV.*tempCOV);
    dacov = dacov.*tempCOV;
    end
else
  decouple = -1;
  end

figure(6); clf; pcolor(log10(abs(double(dacov)))); colorbar; 
title('covariance matr'); pause(0.1)
figure(7); clf; plot(diag(dacov)); pause(0.1); title('diag(cov)')