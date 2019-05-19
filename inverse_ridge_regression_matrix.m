function [invmatr,newmatr] = inverse_ridge_regression_matrix(oldmatr,kmax);

junk = (oldmatr);
[U,S,V] = svd(junk);
diagS = diag(S);
fprintf(1,'max/min eigenvalues before RR = %8.6e %8.6e rank %8.6e \n',max(diagS),min(diagS),max(diagS)/min(diagS));

delta = (max(diagS)-min(diagS)*kmax)/(kmax-1);
[Unew,Snew,Vnew] = svd(junk + delta * eye(size(junk)));  
diagSnew = diag(Snew);
fprintf(1,'max/min eigenvalues after  RR = %8.6e %8.6e rank %8.6e \n',max(diagSnew),min(diagSnew),max(diagSnew)/min(diagSnew));

junk = junk + delta * eye(size(junk));
newmatr = junk;

invmatr = inv(junk);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
