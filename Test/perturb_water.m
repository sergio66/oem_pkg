function p1 = perturb_water(huse,puse,perturb,raPlevs)

% function p1 = perturb_water(puse,perturb,iaLevs)
% takes in rh and applies perturbation as needed; makes sure RH <= 99.99 
%   assumes perturbation in fraction
%   puse : pressures in mb, gas_1 in g/g
%   raPlevs : between which two pressure levels to perturb (-1 = all)

if huse.glist(1) ~= 1
  error('assumes first gas = 1')
end

if huse.gunit(1) ~= 21
  error('assumes units = g/g')
end

if huse.ptype > 0
  error('assumes h.ptype = 0')
end

raPlevs = sort(raPlevs);

p1 = puse;

for ii = 1 : length(puse.stemp)
  levs = 1:puse.nlevs(ii);
  rh = mmr2rh(puse.plevs(levs,ii),puse.ptemp(levs,ii),puse.gas_1(levs,ii));

  rh = rh*(1+perturb);
  oo = find(rh > 100); rh(oo) = 99.99; 
  if length(oo) > 0
    disp('warning : perturbation made a few layers have >100 percent humidity; resetting to 99.999 ')
  end

  %% change RH to ppmv
  ynew = toppmv(puse.plevs(levs,ii),puse.ptemp(levs,ii),rh,18,40);

  % change ppmv back to MMR g/g
  mdair = 28.966;  % molecular mass of dry air
  mass_g = 18;     % molecular mass of WV
  rjunk =1E+6 * mdair/mass_g;
  ynewx = ynew/rjunk;

  if raPlevs == -1
    p1.gas_1(levs,ii) = ynewx;
  else
    woosh = find(puse.plevs(levs,ii) >= raPlevs(1) & puse.plevs(levs,ii) <= raPlevs(2));
    if length(woosh) >= 1
      p1.gas_1(levs(woosh),ii) = ynewx(woosh);
    end
  end

end
