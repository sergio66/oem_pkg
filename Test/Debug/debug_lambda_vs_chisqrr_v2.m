figure(9);
% $$$   scatter(chisqrr,x00diff,60,1:length(lambda),'filled'); colorbar; grid
%% bug in scatter if only have 3 chisqrr's
plot(chisqrr,x00diff,'+-');
title('l curve'); xlabel('|y-kx|^{2}=\chi^2');    ylabel('|x0-x|^{2}');
% $$$   hold on; plot(chisqrr,x00diff); hold off; 
disp('these are the co2 oem retrievals for different lambdas')
coeffsr(:,1)'*2.2
