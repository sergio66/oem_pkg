function [invmatr,newmatr] = inverse_minimum_eigenvalue_matrix_optim(oldmatr,krange,sigrange,comment);

for ii = 1 : length(krange)
  for jj = 1 : length(sigrange)
    [junk1,junk2] = inverse_minimum_eigenvalue_matrix(oldmatr,krange(ii),sigrange(jj),' ');
    junk = (junk1*oldmatr + oldmatr*junk1)/2;
    chisqr(ii,jj) = norm(eye(size(oldmatr)) - junk,'fro');
  end
end

contour(log10(krange),log10(sigrange),log10(chisqr')); colorbar; colormap jet
xlabel('krange = rank'); ylabel('sigrange = min eig');
if nargin == 4
  title(comment);
end

%% now find minimum!!!
[ii,jj] = find(chisqr == min(chisqr));
ii = ii(1);
jj = jj(1);
fprintf(1,'  optimize inverse_minimum_eigenvalue for %s : rank = %8.6e min eig %8.6e \n',comment,krange(ii),sigrange(jj));
[invmatr,newmatr] = inverse_minimum_eigenvalue_matrix(oldmatr,krange(ii),sigrange(jj),comment);

disp(' ')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
