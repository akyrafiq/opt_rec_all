%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%       otimize resource allocation for mm-1 cost-function
%	z == Phi_0 Lagrange parameter controling sum of loads
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
global N M0 la z z0 mui mu
%--------------------------------------------------------------------
%  initialization
%--------------------------------------------------------------------
N  = 10000;         % number of processes
M0 =     1;		    % exponentially distibituted resources MEAN
mui = rand(1,N);	% exprnd(M0,1,N);%  resource parameters
mu  = sum(mui);		% total capacity
%--------------------------------------------------------------------
function y = f (z)
global mu mui la
  y1= mu - la -sum(sqrt(mui/(z*la)).*(mui>sqrt(mui/(z*la))));
  y = y1 -sum(mui.*(mui < sqrt(mui/(z*la))));
end

mdat= [];
x= 0.01:0.01:0.99;
Nl=length(x);
z0=0.2;

for i=1:Nl
  xi=x(i)
  la=x(i)*mu;
  [z, fval, info] = fsolve("f",z0)
  lai = max(mui - sqrt(mui/(z*la)),0);
  D = sum(lai./(mui-lai))/la;
  mdat = [mdat; x(i), z, sum((mui < sqrt(mui/(z*la)))),D];
  z0=z
end

save DRouting/ualloc2  mdat; % output distr of gc probs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%		that's it
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
