function [invmatr,newmatr] = inverse_minimum_eigenvalue_matrix(oldmatr,kmax,sigmin,comment);

warning off

iPrint = -1;
if nargin == 4
  if length(comment) > 1
    iPrint = 1;
  elseif length(comment) == 1 & comment == ' '
    iPrint = -1;
  end
end

junk = (oldmatr);
[U,S,V] = svd(junk);
diagS = diag(S);
%if nargin == 3
%  fprintf(1,'max/min eigenvalues before ME = %8.6e %8.6e rank %8.6e \n',max(diagS),min(diagS),max(diagS)/min(diagS));
if iPrint > 0
  fprintf(1,'max/min eigenvalues before ME = %8.6e %8.6e rank %8.6e for %s \n',max(diagS),min(diagS),max(diagS)/min(diagS),comment);
end

iVers = 1; %% orig code, mu buggy mistake
iVers = 2; %% closer to what is in ImproveConditioning_UKMO.pdf

if iVers == 1
  %% orig code, version 1
  kmax = 1;   %% testing
  small = find(diagS < sigmin);
  diagS(small) = sigmin;
  S = diag(diagS);
  junk = U * S * V';

elseif iVers == 2
  %% newcode, version 2
  lambda1 = max(diagS);
  T1   = lambda1/kmax;
  small = find(diagS < T1);
  diagS(small) = T1;
  diagS(1)     = lambda1;
  S = diag(diagS);
  junk = U * S * V';
end

%if nargin == 3
%  fprintf(1,'max/min eigenvalues after  ME = %8.6e %8.6e rank %8.6e \n',max(diagS),min(diagS),max(diagS)/min(diagS));
if iPrint > 0
  fprintf(1,'max/min eigenvalues after  ME = %8.6e %8.6e rank %8.6e for %s \n',max(diagS),min(diagS),max(diagS)/min(diagS),comment);
end

newmatr = junk;
invmatr = inv(junk);

warning on

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
