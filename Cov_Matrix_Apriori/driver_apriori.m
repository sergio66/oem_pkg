QSTjacindex = [1 1 1 1 0 1];   %% use[CO2 O3 N2O CH4     stemp];
QSTjacindex = [1 1 1 1 1 1];   %% use[CO2 O3 N2O CH4 CFC stemp];

WVjacindex = 999;
if WVjacindex(1) == 999
  WVjacindex = 1:97;
else
  WVjacindex = unique(WVjacindex);
  end

Tjacindex = 999;
if Tjacindex(1) == 999
  Tjacindex = 1:97;
else
  Tjacindex = unique(Tjacindex);
  end

iiBin = input('enter latbin : ');

disp(' ')
  [junk,ERArates,ERAstd,ERArenorm,cutoffs] = ...
              find_cov_damp2(QSTjacindex,WVjacindex,Tjacindex,iiBin,1,-1);

save_apriori.apriori = ERArates;
save_apriori.diagcov = ERAstd;
save_apriori.string  = ...
  'ind indices; from RATES_NEW3/TOOLS/COV_MATRIX_INDIVIDUAL/driver_apriori.m';
save_apriori.QSTjacindex = QSTjacindex;
save_apriori.WVjacindex  = WVjacindex;
save_apriori.Tjacindex   = Tjacindex;
save_apriori.Ind_or_Combined = +1;

saver = ['save ../JUNK_MATRICES/apriori_latbin_' num2str(iiBin) ];
saver = [saver '_individual.mat save_apriori'];
eval(saver);