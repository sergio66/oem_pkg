function [y,jac_new] = run_klayers_sarta(driver,xn,rtpop,rtprp,iIterate,iNewJac);

%% driver = structure which code is dealing with; note it contains "headstruct" and "profstruct" for rtp files
%% rtpop  = .op.rtp file name (temporary) code will use
%% rtprp  = .rp.rtp file name (temporary) code will use
%%    xn  = current atmospheric profile to be put in place of what is in "profstruct"
%% iIterate = integrer, telling which iteration this is 
%% iNewJac  = -1 for no, +1 for yes

jac_new = [];

loader = ['load ' driver.oem.headstruct]; eval(loader); h = h_oem;
loader = ['load ' driver.oem.profstruct]; eval(loader); p = p_oem;

if iIterate == 0
  %% just blindly trust what came in 
  rtpwrite(rtpop,h,[],p,[]);
else
  %% bah humbug ... figure things out, only setup for {STEMP}{WV97}{TEMP97}
  p.stemp = xn(1);
  iOffSet = 1;
  for ii = 1 : driver.jacobian.numlays
    wv = 10^(xn(iOffSet + ii));
    p.gas_1(ii) = wv;
  end
  iOffSet = 1 + driver.jacobian.numlays;
  for ii = 1 : driver.jacobian.numlays
    t = xn(iOffSet + ii);
    p.ptemp(ii) = t;
  end
  ix = 1 : driver.jacobian.numlays;
  figure(4); plot(p.ptemp(ix)- p_oem.ptemp(ix), ix, p.gas_1(ix) ./ p_oem.gas_1(ix),ix,'r'); grid
    mmw0 = mmwater_rtp(h,p_oem);
    mmwn = mmwater_rtp(h,p);
    fprintf(1,'mmw ratio n : orig = %8.6f \n',mmwn/mmw0);

    set(gca,'ydir','reverse'); hl = legend('temp(n)-temp(0)','WV(n)/WV(0)'); set(hl,'fontsize',10)

  figure(1); plot(p.ptemp(ix),ix,'r',p_oem.ptemp(ix),ix,'b'); set(gca,'ydir','reverse')
  figure(3); semilogx(p.gas_1(ix),ix,'r',p_oem.gas_1(ix),ix,'b'); set(gca,'ydir','reverse')
  rtpwrite(rtpop,h,[],p,[]);
  figure(2);
end

sarta = driver.oem.sarta;
sartaer   = ['!' sarta   ' fin=' rtpop '  fout=' rtprp]; eval(sartaer);
[hh,hha,pp,ppa] = rtpread(rtprp);
y = rad2bt(hh.vchan,pp.rcalc);

if iNewJac > 0

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