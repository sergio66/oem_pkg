%%
%% called by perturb_49regr_profiles.m to save stuff for retrievals
%%
%% see /strowdataN/home/sergio/MATLABCODE/oem_pkg_run/Generic/override_defaults_generic.m
%% what you want to retrieve should be in a matlab which has variables of following form
%%
%%  b_bias          36x2378x10            6848640  double
%%  b_cal           36x2378x10            6848640  double
%%  b_err_bias      36x2378x10            6848640  double
%%  b_err_cal       36x2378x10            6848640  double
%%  b_err_obs       36x2378x10            6848640  double
%%  b_obs           36x2378x10            6848640  double
%%
%% where what we want to retrieve are in the SECOND field ie b_obs(:,:.2)
%% where what we want to retrieve are in the SECOND field ie b_obs(:,:.2)
%% where what we want to retrieve are in the SECOND field ie b_obs(:,:.2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%>>>>>>>>>>>>> for set_up_oem_pkg_spectraljacs1
pxPERT(1) = h2p(3000);
pxPERT(2) = 1200;
p49new = perturb_water(h49,p49,+0.10,pxPERT);
rtpwrite('junk.ip.rtp',h49,ha49,p49new,pa49);
klayerser = ['!' klayers ' fin=junk.ip.rtp fout=junk.op.rtp']; eval(klayerser);
sartaer   = ['!' sarta   ' fin=junk.op.rtp fout=junk.rp.rtp']; eval(sartaer);
[h49x,ha49x,p49newy,pa49x] = rtpread('junk.rp.rtp');
%%%%%%%%%%%%%%%%%%%%%%%%%>>>>>>>>>>>>>

plot(h49x.vchan,rad2bt(h49x.vchan,p49x.rcalc))
plot(h49x.vchan,rad2bt(h49x.vchan,p49x.rcalc)-rad2bt(h49x.vchan,p49newx.rcalc))
plot(h49x.vchan,rad2bt(h49x.vchan,p49x.rcalc)-rad2bt(h49x.vchan,p49newy.rcalc))

iaCold = find(p49.stemp <= 273);
iaMed  = find(p49.stemp > 273 & p49.stemp <= 290);
iaHot  = find(p49.stemp > 290); 

plot(h49x.vchan,rad2bt(h49x.vchan,p49x.rcalc(:,iaCold))-rad2bt(h49x.vchan,p49newx.rcalc(:,iaCold)),'b',...
     h49x.vchan,rad2bt(h49x.vchan,p49x.rcalc(:,iaMed))-rad2bt(h49x.vchan,p49newx.rcalc(:,iaMed)),'g',...
     h49x.vchan,rad2bt(h49x.vchan,p49x.rcalc(:,iaHot))-rad2bt(h49x.vchan,p49newx.rcalc(:,iaHot)),'r')
xlabel('Wavenumber cm-1'); ylabel('\delta BT(K)');
title('cold(b) warm (g) and hot (r) profiles')

t49x    = rad2bt(h49x.vchan,p49x.rcalc);
t49newx = rad2bt(h49x.vchan,p49newx.rcalc);
t49newy = rad2bt(h49x.vchan,p49newy.rcalc);
%{
plot(h49x.vchan,t49x(:,iaCold)-t49newy(:,iaCold),'b',...
     h49x.vchan,t49x(:,iaMed)-t49newy(:,iaMed),'g',...
     h49x.vchan,t49x(:,iaHot)-t49newy(:,iaHot),'r',...
     h.vchan(chans41),allspectra0 - allspectra10,'k')
plot(h49x.vchan,nanmean(t49x(:,iaCold)'-t49newy(:,iaCold)'),'b',...
     h49x.vchan,nanmean(t49x(:,iaMed)'-t49newy(:,iaMed)'),'g',...
     h49x.vchan,nanmean(t49x(:,iaHot)'-t49newy(:,iaHot)'),'r',...
     h.vchan(chans41),allspectra0 - allspectra10,'k')
%}
hl = legend('cold','warm','hot','rico'); grid
set(hl,'fontsize',10); grid
title('Prof0 - +0.1Pert from 0.3 km'); xlabel('Wavenumber cm-1'); ylabel('\delta BT(K)')
axis([650 2700 -0.001 +0.4])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% these are the spectral observations, with the perturbations to WV in lowest 1.5 km
nedt = instr_chans('airs',2);
clear b_*
for ii = 1 : 1
  b_obs(ii,:,2) = t49newy(:,ii)';
  b_err_obs(ii,:,2) = nedt;
end
comment = 'see /home/sergio/MATLABCODE/LES/perturb_49regr_profiles.m';
save /strowdataN/home/sergio/MATLABCODE/oem_pkg_run/Generic/test_prof.mat b_obs b_err_obs comment

clear apriori*
apriori(1)       = p49newx.stemp(1) + 1.00;
apriori(002:098) = log10(p49newx.gas_1(1:97,1)*1.01);
apriori(099:195) = p49newx.ptemp(1:97,1) + 1.00;
apriori = apriori';
save /strowdataN/home/sergio/MATLABCODE/oem_pkg_run/Generic/test_prof_apriori.mat apriori comment

[h_oem,p_oem] = replicate_rtp_headprof(h49x,p49x,1,1);
p_oem.stemp = p_oem.stemp + 1;
p_oem.ptemp = p_oem.ptemp + 1;
p_oem.gas_1 = p_oem.gas_1 * 1.01;
save /strowdataN/home/sergio/MATLABCODE/oem_pkg_run/Generic/test_prof_headstruct.mat h_oem
save /strowdataN/home/sergio/MATLABCODE/oem_pkg_run/Generic/test_prof_profstruct.mat p_oem

%%%%%%%%%%%%%%%%%%%%%%%%%

axS = load('/asl/s1/sergio/AIRSPRODUCTS_JACOBIANS/STD/surface_jac.mat');
axW = load('/asl/s1/sergio/AIRSPRODUCTS_JACOBIANS/TRP/g1_jac.mat');
axT = load('/asl/s1/sergio/AIRSPRODUCTS_JACOBIANS/TRP/temp_jac.mat');

BST = ones(1,1)*(1);             %% error in ECMWF T(z) about 1.0
BQ = ones(1,97)*log10(1 + 0.1);  %% error in ECMWF q(z) about 0.1
BT = ones(1,97)*(1);             %% error in ECMWF T(z) about 1.0
B = [BST BT BQ]; B = diag(B); B = B.*B;  

clear B
hgts = p2h(p_oem.plevs);
hgts(1) = hgts(2) + abs(hgts(2)-hgts(3))*1.5;
B = zeros(195,195);
B(1,1) = 1^2;
clear bonk
%% also tie the surface to the uper levels
scalehgt = 100;
for ii = 1 : 97
  jj = 97;
  junk = abs(hgts(ii)-hgts(jj))/scalehgt;
  bonk(ii) = exp(-junk);
  junk = abs(ii-jj)/10;
  bonkij(ii) = exp(-junk);
end
B(1,99:195) = bonk;
B(99:195,1) = bonk';

clear bonk
scalehgt = 5000;
scalehgt = 2500;
scalehgt = 1000;
for ii = 1 : 97
  for jj = 1 : 97
    junk = abs(hgts(ii)-hgts(jj))/scalehgt;
    bonk(ii,jj) = exp(-junk);
  end
end
B(002:098,002:098) = (log10(1 + 0.1))^2 * bonk;
B(099:195,099:195) = (1.0^2) * bonk;

cov = ones(1,195);
cov = diag(cov);
cov = (B);
save /strowdataN/home/sergio/MATLABCODE/oem_pkg_run/Generic/test_prof_cov.mat cov

%%%%%%%%%%%%%%%%%%%%%%%%%

if isfield(p49x,'landtype')
 p49x = rmfield(p49x,'landtype');
end
if isfield(p49x,'freqcal')
  p49x = rmfield(p49x,'freqcal');
end
if isfield(p49x,'robsqual')
  p49x = rmfield(p49x,'robsqual');
end
if isfield(p49x,'ifov')
  p49x = rmfield(p49x,'ifov');
end
p49x.nemis = single(p49x.nemis);
p49x.nlevs = single(p49x.nlevs);
p49x.ctype2 = single(p49x.ctype2);
p49x.upwell = single(p49x.upwell);
p49x.clrflag = single(p49x.clrflag);
p49x.ctype = single(p49x.ctype);
p49x.findex = single(p49x.findex);
p49x.atrack = single(p49x.atrack);
p49x.xtrack = single(p49x.xtrack);
%p49x.ifov = single(p49x.ifov);
p49x.pnote = single(p49x.pnote);
p49x.iudef = single(p49x.iudef);
p49x.itype = single(p49x.itype);
%p49x.robsqual = single(p49x.robsqual);

[h1,p1] = replicate_rtp_headprof(h49x,p49x,1,97*2+1+1);
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
rtpwrite('junk.op.rtp',h1,ha49,p1,pa49);
%klayerser = ['!' klayers ' fin=junk.ip.rtp fout=junk.op.rtp']; eval(klayerser);
sartaer   = ['!' sarta   ' fin=junk.op.rtp fout=junk.rp.rtp']; eval(sartaer);
[h1x,ha49x,p1x,pa49x] = rtpread('junk.rp.rtp');

tpert = rad2bt(h1x.vchan,p1x.rcalc);
tpert = tpert(:,2:196) - tpert(:,1)*ones(1,195);
jj = 1; %% stemp
  M_TS_jac_all(1,:,1)       = tpert(:,jj)/1;
for ii = 1 : 97  %% wv97
  jj = jj + 1;
  %% [f(q(1+delta) - f(q)]/[log10(q(1+delta)) - log10(q)] = 
  %% [f(q(1+delta) - f(q)]/[log10(q) + log10(1+delta) - log10(q)] = 
  %% [f(q(1+delta) - f(q)]/[log10(1+delta)] = 
  M_TS_jac_all(1,:,jj) = tpert(:,jj)/log10(1.1);
end
for ii = 1 : 97  %% ptemp 97
  jj = jj + 1;
  M_TS_jac_all(1,:,ii+97+1) = tpert(:,jj)/1;
end
qrenorm = ones(1,97*2+1);
f = axT.fout;
save /strowdataN/home/sergio/MATLABCODE/oem_pkg_run/Generic/test_prof_jac.mat M_TS_jac_all f comment qrenorm


