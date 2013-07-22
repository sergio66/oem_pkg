function [y,jac_new] = run_klayers_sarta(driver,xn,rtpop,rtprp,iIterate,iNewJac);

%% driver = structure which code is dealing with; note it contains "headstruct" and "profstruct" for rtp files
%% rtpop  = .op.rtp file name (temporary) code will use
%% rtprp  = .rp.rtp file name (temporary) code will use
%%    xn  = current atmospheric profile to be put in place of what is in "profstruct"
%% iIterate = integrer, telling which iteration this is 
%% iNewJac  = -1 for no, +1 for yes

addpath /asl/matlib/rtptools/
addpath /asl/matlib/h4tools/

jac_new = [];

loader = ['load ' driver.oem.headstruct]; eval(loader); h = h_oem;
loader = ['load ' driver.oem.profstruct]; eval(loader); p = p_oem;

[h,p_oem] = subset_rtp(h,p_oem,[],[],driver.iibin);
[h,p]     = subset_rtp(h,p,    [],[],driver.iibin);

if iIterate <= 0
  %% just blindly trust what came in 
  rtpwrite(rtpop,h,[],p,[]);
else
  %% bah humbug ... figure things out

  %% this first bit does SARTA variable tracegases, stemp and scanang, and cld vars
  str_allowedU = {'CO2' 'O3' 'N2O' 'CO' 'CH4' 'SO2' 'HNO3' 'STEMP' 'SCANANG' 'CPSIZE1' 'CNGWAT1' 'CPRTOP1' 'CFRAC1' 'CPSIZE2' 'CNGWAT2' 'CPRTOP2' 'CFRAC2'};
  str_allowedL = {'co2' 'o3' 'n2o' 'co' 'ch4' 'so2' 'hno3' 'stemp' 'scanang' 'cpsize1' 'cngwat1' 'cprtop1' 'cfrac1' 'cpsize2' 'cngwat2' 'cprtop2' 'cfrac2'};
  for gg = 1 : length(driver.jacobian.qstnames)
    iX = -1;
    dastr = driver.jacobian.qstnames{gg};
    iCnt = 0;
    while iX < 0 & iCnt < length(str_allowedU)
      iCnt = iCnt + 1;
      xstrU = str_allowedU{iCnt};
      xstrL = str_allowedL{iCnt};
      if strcmp(dastr,xstrU) | strcmp(dastr,xstrL)
        iX = +1;
        if iCnt == 1
          gX = 2;    %% column CO2
          p.gas_2 = p.gas_2 * 10^(xn(gg));
          fprintf(1,'CO2 mult at iteration # %2i = %8.6f \n',iIterate,10^(xn(gg)))
        elseif iCnt == 2
          gX = 3;    %% column O3
          p.gas_3 = p.gas_3 * 10^(xn(gg));
          fprintf(1,'O3  mult at iteration # %2i = %8.6f \n',iIterate,10^(xn(gg)))
        elseif iCnt == 3
          gX = 4;    %% column N2O
          p.gas_4 = p.gas_4 * 10^(xn(gg));
          fprintf(1,'N2O mult at iteration # %2i = %8.6f \n',iIterate,10^(xn(gg)))
        elseif iCnt == 4
          gX = 5;    %% column CO
          p.gas_5 = p.gas_5 * 10^(xn(gg));
          fprintf(1,'CO  mult at iteration # %2i = %8.6f \n',iIterate,10^(xn(gg)))
        elseif iCnt == 5
          gX = 6;    %% column CH4
          p.gas_6 = p.gas_6 * 10^(xn(gg));
          fprintf(1,'CH4 mult at iteration # %2i = %8.6f \n',iIterate,10^(xn(gg)))
        elseif iCnt == 6
          gX = 9;    %% column SO2
          p.gas_9 = p.gas_9 * 10^(xn(gg));
          fprintf(1,'SO2 mult at iteration # %2i = %8.6f \n',iIterate,10^(xn(gg)))
        elseif iCnt == 7
          gX = 12;   %% column HNO3
          p.gas_12 = p.gas_12 * 10^(xn(gg));
          fprintf(1,'HN03 mult at iteration # %2i = %8.6f \n',iIterate,10^(xn(gg)))
        elseif iCnt == 8
          gX = -1;   %% stemp
          p.stemp = xn(gg);
          fprintf(1,'stemp offset at iteration # %2i = %8.6f, final = %8.6f \n',iIterate,xn(gg)-p_oem.stemp,p.stemp)
        elseif iCnt == 9
          gX = -2;   %% scanang
          p.scanang = xn(gg);
          p.satzen  = xn(gg);
          fprintf(1,'scanang offset at iteration # %2i = %8.6f, final = %8.6f \n',iIterate,xn(gg)-p_oem.scanang,p.scanang)
        elseif iCnt == 10
          gX = -3;   %% cpsize1
          p.cpsize = p.cpsize * 10^(xn(gg));
          fprintf(1,'cpsize1 mult at iteration # %2i = %8.6f, final = %8.6f \n',iIterate,10^(xn(gg)),p.cpsize)
        elseif iCnt == 11
          gX = -4;   %% cngwat1
          p.cngwat = p.cngwat * 10^(xn(gg));
          fprintf(1,'cngwat1 mult at iteration # %2i = %8.6f, final = %8.6f \n',iIterate,10^(xn(gg)),p.cngwat)
        elseif iCnt == 12
          gX = -5;   %% cprtop1
          p.cprtop = min(p.cprtop * 10^(xn(gg)),p.cprbot-10);
          fprintf(1,'cprtop1 mult at iteration # %2i = %8.6f, final = %8.6f \n',iIterate,10^(xn(gg)),p.cprtop)
        elseif iCnt == 13
          gX = -6;   %% cfrac1
          p.cfrac = max(p.cfrac * 10^(xn(gg)),1);
          fprintf(1,'cfrac1 mult at iteration # %2i = %8.6f, final = %8.6f \n',iIterate,10^(xn(gg)),p.cfrac)
        elseif iCnt == 14
          gX = -3;   %% cpsize2
          p.cpsize2 = p.cpsize2 * 10^(xn(gg));
          fprintf(1,'cpsize2 mult at iteration # %2i = %8.6f, final = %8.6f \n',iIterate,10^(xn(gg)),p.cpsize2)
        elseif iCnt == 15
          gX = -4;   %% cngwat2
          p.cngwat2 = p.cngwat2 * 10^(xn(gg));
          fprintf(1,'cngwat2 mult at iteration # %2i = %8.6f, final = %8.6f \n',iIterate,10^(xn(gg)),p.cngwat2)
        elseif iCnt == 16
          gX = -5;   %% cprtop2
          p.cprtop2 = min(p.cprtop2 * 10^(xn(gg)),p.cprbot2-10);
          fprintf(1,'cprtop2 mult at iteration # %2i = %8.6f, final = %8.6f \n',iIterate,10^(xn(gg)), p.cprtop2)
        elseif iCnt == 17
          gX = -6;   %% cfrac2
          p.cfrac2 = max(p.cfrac2 * 10^(xn(gg)),1);
          fprintf(1,'cfrac2 mult at iteration # %2i = %8.6f, final = %8.6f \n',iIterate,10^(xn(gg)), p.cfrac2)
        end   %% if iCnt
      end     %% if strcmp 
    end       %% while iX < 0
    if iX < 0
      fprintf(1,'\n >>>> driver.jacobian.qstnames %3i string %s \n',gg,dastr)
      disp('>>>>> warning could not be matched to SARTA profile fields, not a standard product code deals with')
    end
  end         %% for gg
  %% do gases
  for gg = 1 : driver.jacobian.numQlays
    gasID = driver.jacobian.gasID_Qlays(gg);
    fprintf(1,'layers at iteration # %2i for gasID %2i \n',iIterate,gasID)
    %% which gases do we want to retrieve
    iOffSet = length(driver.jacobian.qstYesOrNo) + (gg-1)*driver.jacobian.numlays;
    for ii = 1 : driver.jacobian.numlays
      gas = 10^(xn(iOffSet + ii));
      str = ['p.gas_' num2str(gasID) '(ii) = gas;'];
      eval(str)
    end
  end

  %% finally do temperature
  fprintf(1,'layers at iteration # %2i for T(z) \n',iIterate)
  iOffSet = length(driver.jacobian.qstYesOrNo) + driver.jacobian.numlays * driver.jacobian.numQlays;
  for ii = 1 : driver.jacobian.numlays
    t = xn(iOffSet + ii);
    p.ptemp(ii) = t;
  end
  ix = 1 : driver.jacobian.numlays;

  figure(1); plot(p.ptemp(ix),ix,'r',p_oem.ptemp(ix),ix,'b'); set(gca,'ydir','reverse')
  figure(3); semilogx(p.gas_1(ix),ix,'r',p_oem.gas_1(ix),ix,'b'); set(gca,'ydir','reverse')
  figure(4); plot(p.ptemp(ix)- p_oem.ptemp(ix), ix, p.gas_1(ix) ./ p_oem.gas_1(ix),ix,'r'); grid
    mmw0 = mmwater_rtp(h,p_oem);
    mmwn = mmwater_rtp(h,p);
    fprintf(1,'mmw ratio n : orig = %8.6f \n',mmwn/mmw0);

    set(gca,'ydir','reverse'); hl = legend('temp(n)-temp(0)','WV(n)/WV(0)'); set(hl,'fontsize',10)

  rtpwrite(rtpop,h,[],p,[]);
  figure(2); pause(0.1)
end

sarta = driver.oem.sarta;
sartaer   = ['!' sarta   ' fin=' rtpop '  fout=' rtprp]; eval(sartaer);
[hh,hha,pp,ppa] = rtpread(rtprp);
y = rad2bt(hh.vchan,pp.rcalc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if iNewJac > 0

 if driver.jacobian.numQlays > 0
   error('arrrr only have this coded up for WV and temperature profiles')
 end

 if isfield(p,'landtype')
   p = rmfield(p,'landtype');
  end
  if isfield(p,'freqcal')
    p = rmfield(p,'freqcal');
  end
  if isfield(p,'robsqual')
    p = rmfield(p,'robsqual');
  end
  if isfield(p,'ifov')
    p = rmfield(p,'ifov');
  end
  %{
  p.nemis = single(p.nemis);
  p.nlevs = single(p.nlevs);
  p.ctype2 = single(p.ctype2);
  p.upwell = single(p.upwell);
  p.clrflag = single(p.clrflag);
  p.ctype = single(p.ctype);
  p.findex = single(p.findex);
  p.atrack = single(p.atrack);
  p.xtrack = single(p.xtrack);
  %p.ifov = single(p.ifov);
  p.pnote = single(p.pnote);
  p.iudef = single(p.iudef);
  p.itype = single(p.itype);
  %p.robsqual = single(p.robsqual);
  %}

  [h1,p1] = replicate_rtp_headprof(h,p,1,97*2+1+1);
  jj = 1; %% do not change, the control
  jj = 2; 
    p1.stemp(jj) = p1.stemp(jj)+1; %% stemp;
  for ii = 1 : 97
    jj = jj + 1;
    p1.gas_1(ii,jj) =   p1.gas_1(ii,jj) * 1.1;
  end
  for ii = 1 : 97
    jj = jj + 1;
    p1.ptemp(ii,jj) =   p1.ptemp(ii,jj) + 1;
  end
  rtpwrite('junk.op.rtp',h1,[],p1,[]);
  %klayerser = ['!' klayers ' fin=junk.ip.rtp fout=junk.op.rtp']; eval(klayerser);
  sartaer   = ['!' sarta   ' fin=junk.op.rtp fout=junk.rp.rtp']; eval(sartaer);
  [h1x,ha49x,p1x,pa49x] = rtpread('junk.rp.rtp');

  tpert = rad2bt(h1x.vchan,p1x.rcalc);
  tpert = tpert(:,2:196) - tpert(:,1)*ones(1,195);
  jj = 1; %% stemp
    jac_new(:,1)       = tpert(:,jj)/1;
  for ii = 1 : 97  %% wv97
    jj = jj + 1;
    %% [f(q(1+delta) - f(q)]/[log10(q(1+delta)) - log10(q)] = 
    %% [f(q(1+delta) - f(q)]/[log10(q) + log10(1+delta) - log10(q)] = 
    %% [f(q(1+delta) - f(q)]/[log10(1+delta)] = 
    jac_new(:,jj) = tpert(:,jj)/log10(1.1);
  end
  for ii = 1 : 97  %% ptemp 97
    jj = jj + 1;
    jac_new(:,ii+97+1) = tpert(:,jj)/1;
  end
end

rmer = ['!/bin/rm ' rtpop ' ' rtprp];
eval(rmer);