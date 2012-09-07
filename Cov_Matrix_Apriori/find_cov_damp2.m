function [dacov,ecmwf_rates,ecmwf_rates_unc,renorm,cutoffs] = ...
  find_cov_damp2(driver,iDoCov)

% grows a exp(-((i-j)/damp)^2) for W,T and diag for Q
% same as find_cov_damp except it loads in ECMWF profiles and puts in the
% rates from them

loader = ['load ' driver.ecmfile]; eval(loader)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% trace gases
xozonerate    = sum(ecmwf.ozonerate')/97;  
xstdozonerate = sum(ecmwf.ozoneratestd')/97;  
xozonerate    = max(ecmwf.ozonerate')/97;  
xstdozonerate = max(ecmwf.ozoneratestd')/97;  

%% [CO2  O3   N2O  CH4 CFC]
%% [ppmv frac ppbv ppb ppb] yr-1
ecmwf_rates_gas     = [2.2/370 xozonerate(latBin)    1.0/300 5/1800 -1/1300]; 
ecmwf_rates_unc_gas = [0.2/370 xstdozonerate(latBin) 0.2/300 2/1800 0.2/1300]; 
renormQ = driver.qrenorm(1:6);  %[trace gases and stemp]
ecmwf_rates_gas     = ...
  [2.2 xozonerate(latBin)    1.0 5 -1.0 ecmwf.stemprate(latBin)];
ecmwf_rates_unc_gas = ...
  [0.2 xstdozonerate(latBin) 0.2 2 0.2 ecmwf.stempratestd(latBin)]; 

if sum(QSTjacindex) >= 6
  ix = ones(1,6);
  qCOV = diag(ix);
  ecmwf_rates     = ecmwf_rates_gas./renormQ;
  ecmwf_rates_unc = ecmwf_rates_unc_gas./renormQ;
  qCOV = diag(ecmwf_rates_unc.^2);
else
  ix = ones(1,sum(QSTjacindex));
  qCOV = diag(ix);
  ecmwf_rates     = ecmwf_rates_gas(QSTjacindex)./renormQ(QSTjacindex);
  ecmwf_rates_unc = ecmwf_rates_unc_gas(QSTjacindex)./renormQ(QSTjacindex);
  qCOV = diag(ecmwf_rates_unc.^2);
  end

renorm = renormQ;
cutoffs(1) = length(renorm);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% wv

clear blah sblah
xwaterrate    = ecmwf.waterrate(latBin,:);
xstdwaterrate = ecmwf.waterratestd(latBin,:);
for ix = 1 : length(WVjacindex)
  woof = WVjacindex(ix);
  blah(ix)  = xwaterrate(woof);
  sblah(ix) = xstdwaterrate(woof);
end
xwaterrate    = blah;
xstdwaterrate = sblah;

xwaterrate    = xwaterrate/0.01;
xstdwaterrate = xstdwaterrate/0.01;

renorm = [renorm 0.01*ones(size(xwaterrate))];
cutoffs(2) = length(renorm);

ix = xstdwaterrate;
ixA = ones(length(WVjacindex),1) * ix;
ixB = ixA';
gas1COV = ixA.*ixB;

%% put damps
ix = 1 : length(WVjacindex);
ixA = ones(length(WVjacindex),1) * ix;
ixB = ixA';
xCOV = ixA-ixB;
xCOV = xCOV/covDamp; 
xCOV = exp(-xCOV.*xCOV);

gas1COV = gas1COV.*xCOV; 

ecmwf_rates     = [ecmwf_rates     xwaterrate];
ecmwf_rates_unc = [ecmwf_rates_unc xstdwaterrate];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% temp

clear blah sblah
xtemprate    = ecmwf.ptemprate(latBin,:);
xstdtemprate = ecmwf.ptempratestd(latBin,:);
for ix = 1 : length(Tjacindex)
  woof = Tjacindex(ix);
  blah(ix)  = xtemprate(woof);
  sblah(ix) = xstdtemprate(woof);
  end
xtemprate    = blah;
xstdtemprate = sblah;

xtemprate    = xtemprate/0.1;
xstdtemprate = xstdtemprate/0.1;

renorm = [renorm 0.1*ones(size(xtemprate))];
cutoffs(3) = length(renorm);

ix = xstdtemprate;
ixA = ones(length(xstdtemprate),1) * ix;
ixB = ixA';
tempCOV = ixA.*ixB;

ix = 1 : length(Tjacindex);
ixA = ones(length(Tjacindex),1) * ix;
ixB = ixA';
xCOV = ixA-ixB;
xCOV = xCOV/covDamp; 
xCOV = exp(-xCOV.*xCOV);

tempCOV = tempCOV.*xCOV; 
ecmwf_rates     = [ecmwf_rates     xtemprate];
ecmwf_rates_unc = [ecmwf_rates_unc xstdtemprate];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lala = length(qCOV) + length(gas1COV) + length(tempCOV);
dacov = zeros(lala,lala);
ix = 1:length(qCOV); dacov(ix,ix) = qCOV;
iy = 1:length(gas1COV); ix = length(qCOV)+iy;                 
  dacov(ix,ix) = gas1COV;
iy = 1:length(tempCOV); ix = length(qCOV)+length(gas1COV)+iy; 
  dacov(ix,ix) = tempCOV;

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

ecmwf_rates     = ecmwf_rates';
ecmwf_rates_unc = ecmwf_rates_unc';

if iDoCov > 0
  figure(6); clf; pcolor(dacov); colorbar; 
    title(['damping = ' num2str(covDamp)]); pause(0.1)
  figure(7); clf; plot(diag(dacov)); pause(0.1); title('diag(cov)')
else
  figure(6); clf; plot(ecmwf_rates); 
    title(['apriori rates for latbin ' num2str(latBin)]);
  figure(7); clf; plot(ecmwf_rates_unc); 
    title(['apriori unc in rates for latbin ' num2str(latBin)]);
   pause(0.1);
  end

