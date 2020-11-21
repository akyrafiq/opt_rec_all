%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%   otimize resource allocation for quadratic cost-function
%	z == Phi_0 Lagrange parameter controling sum of loads
%	Determine lambda as function of z =lammbda*phi0 !!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
global N M0 b c la z z0  mui mu
%--------------------------------------------------------------------
%  initialization
%--------------------------------------------------------------------
N  = 1000;              % number of processes
M0 =    1;		% exponentially distibituted resources MEAN
a  =    1;		% parameters of cost function
b  =    1;
c  =    1;
mui = rand(1,N);		% exprnd(M0,1,N);	% resource parameters
mu  = sum(mui);		% total capacity
%--------------------------------------------------------------------
function la = f (z)
global b c mu mui
  x0 = max(4*c*(mui*z-1),0);
  x0 = min((sqrt(b^2 +x0) -b)/(2*c),1);
  la = sum(mui.*x0);
end

mdat= [];
zi= 1:0.1:50;
Nl=length(zi);


for i=1:Nl
  z =zi(i)
  la = f(z)
  if ((la > 0) & (la < mu))
    phi0 = z/la;
    x0i = max(4*c*(z*mui-1),0);
    x0i = min((sqrt(b^2 +x0i) -b)/(2*c),1);			
    D = sum(x0i.*(1+ 0.5*b*x0i +c*x0i.^2/3))/la;
    mdat = [mdat; la/mu, phi0, sum((x0i <=0 )),sum((x0i >= 1)),D];
  end
end

save DRouting/ualloc-111  mdat; % output distr of gc probs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%		that's it
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
