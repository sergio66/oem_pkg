  figure(8); 
  plot(1:length(erastd),erastd.^2,'b*-',1:length(erastd),diag(r0),'r'); 
  title('(b) diag(era uncertainties) (r) diag(cov) ')
  disp('see fig for diag(cov) ...');
  drawnow
