function [y] = run_klayers_sarta(driver,xn,rtpop,rtprp,iIterate);

loader = ['load ' driver.oem.headstruct]; eval(loader); h = h_oem;
loader = ['load ' driver.oem.profstruct]; eval(loader); p = p_oem;

if iIterate == 0
  !! just blindly trust what came in 
  rtpwrite(rtpop,h,[],p,[]);
else
  %% bah humbug ... figure things out, only setup for {STEMP}{WV97}{TEMP97}
  p_oem.stemp = p_oem.stemp + xn(1);
  for ii = 1 : driver.jacobian.numlays
    wv = p_oem.gas_1(ii);
    dwv = exp(xn(1+ii));
    p_oem.gas_1(ii) = wv * dwv;
  end
  for ii = 1 : driver.jacobian.numlays
    t = p_oem.ptemp(ii);
    dt = xn(1+97+ii);
    p_oem.ptemp(ii) = t + dt;
  end
  figure(1); plot(p.ptemp-p_oem.ptemp,1:101, p.gas_1 ./ p_oem.gas_1,1:101)
  p = p_oem;
  figure(2);
end

sarta = driver.oem.sarta;
sartaer   = ['!' sarta   ' fin=' rtpop '  fout=' rtprp]; eval(sartaer);
[hh,hha,pp,ppa] = rtpread(rtprp);
y = rad2bt(hh.vchan,pp.rcalc);

