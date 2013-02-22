function [xn,idx] = despike2(x,fac)

%% see http://marymount.mmm.edu/faculty/kjordahl/matlab/despike.m

%DESPIKE  5-point window data despiking routine
%
%  XN = DESPIKE(X)
%  [XN,IDX] = DESPIKE(X,FAC)
%
%  XN = despiked data
%  IDX = indices of removed spikes
%  FAC = tolerance factor (default 5)
%
% Despike function with 5 point window
% Algorithm of Trehu & Sutton, Mar. Geophys. Res., 16, 91-103, 1994.
%

% Kelsey Jordahl, MBARI
% Time-stamp: <Fri Oct 19 10:50:21 PDT 2001>

if nargin < 2
  % set default tolerance factor
  fac=5; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if mod(fac,2) == 0
  %% make sure fac is odd
  fac = fac + 1;
end

if fac < 5
  fac = 5;
end
  	       
facx = (fac-1)/2;

j = 3:length(x)-2;                        % don't use points < 2 samp. from edge
j = facx+1 :length(x)-facx;               % don't use points < (N-1)/2 samp. from edge
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xn = x;

adiff1 = (abs(x(j-2)-x(j-1)) + abs(x(j+1)-x(j+2)))/2;
adiff2 = (abs(x(j-1)-x(j)) + abs(x(j)-x(j+1)))/2;
idx    = find(adiff2>(fac*adiff1))+2;     % indices of spikes (correct for edges)
xn(idx)=(x(idx+1)+x(idx-1))/2;            % interpolate bad points

% xn(idx)=NaN;                            % throw out bad points
