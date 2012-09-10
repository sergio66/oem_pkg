 
  ii = 0;
  if qstjacindex(1)
    jacobian.co2 = true;
    ii = ii + 1;
    co2(iibin) = oem.finalrates(ii);     dco2(iibin) = oem.finalsigs(ii); 
  else
    jacobian.co2 = false;
    co2(iibin) = NaN;  dco2(iibin) = NaN;
  end

  if qstjacindex(2)
    ii = ii + 1;
    jacobian.co2strat = true;
    co2strat(iibin) = oem.finalrates(ii);     dco2strat(iibin) = oem.finalsigs(ii); 
  else
    jacobian.co2strat = false;
    co2strat(iibin) = NaN;  dco2(iibin) = NaN;
  end

  if qstjacindex(3)
    ii = ii + 1;
    jacobian.o3 = true;
    o3(iibin)  = oem.finalrates(ii);     do3(iibin)  = oem.finalsigs(ii); 
  else
    jacobian.o3 = false;
    o3(iibin) = NaN;  do3(iibin) = NaN;
  end

  if qstjacindex(4)
    ii = ii + 1;
    jacobian.o3strat = true;
    o3strat(iibin)  = oem.finalrates(ii);     do3strat(iibin)  = oem.finalsigs(ii); 
  else
    jacobian.o3strat = false;
    o3strat(iibin) = NaN;  do3strat(iibin) = NaN;
  end

  if qstjacindex(5)
    ii = ii + 1;
    jacobian.n2o = true;
    n2o(iibin) = oem.finalrates(ii);     dn2o(iibin) = oem.finalsigs(ii); 
  else
    jacobian.n2o = false;
    n2o(iibin) = NaN;  dn2o(iibin) = NaN;
  end

  if qstjacindex(6)
    ii = ii + 1;
    jacobian.co = true;
    co(iibin) = oem.finalrates(ii);     dco(iibin) = oem.finalsigs(ii); 
  else
    jacobian.co = false;
    co(iibin) = NaN;  dco(iibin) = NaN;
  end

  if qstjacindex(7)
    ii = ii + 1;
    jacobian.ch4 = true;
    ch4(iibin) = oem.finalrates(ii);     dch4(iibin) = oem.finalsigs(ii); 
  else
    jacobian.ch4 = false;
    ch4(iibin) = NaN;  dch4(iibin) = NaN;
  end

  if qstjacindex(8)
    ii = ii + 1;
    jacobian.cfc = true;
    cfc(iibin) = oem.finalrates(ii);     dcfc(iibin) = oem.finalsigs(ii); 
  else
    jacobian.cfc = false;
    cfc(iibin) = NaN;  dcfc(iibin) = NaN;
  end

  if qstjacindex(9)
    ii = ii + 1;
    jacobian.hdo = true;
    hdo(iibin) = oem.finalrates(ii);     dhdo(iibin) = oem.finalsigs(ii); 
  else
    jacobian.hdo = false;
    hdo(iibin) = NaN;  dhdo(iibin) = NaN;
  end

  if qstjacindex(10)
    ii = ii + 1;
    jacobian.stemp = true;
    stemp(iibin) = oem.finalrates(ii);   dstemp(iibin) = oem.finalsigs(ii); 
  else
    jacobian.stemp = false;
    stemp(iibin) = NaN;  dstemp(iibin) = NaN;
    end
