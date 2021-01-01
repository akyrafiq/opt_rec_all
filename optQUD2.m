%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%  otimize resource allocation for quadratic cost-function
%	z == Phi_0 Lagrange parameter controling sum of loads x(i)==f_i
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
global la z z0 mua na rha mu aa ba
%--------------------------------------------------------------------
%  initialization
%--------------------------------------------------------------------
na = [25, 10, 5];  % number of processes in three classes
mua = [1, 5, 25];  % exponentially distibituted resources MEAN
mu  = sum(na.*mua);% total capacity
aa=1; ba=0.2;
%--------------------------------------------------------------------

%--------------------------------------------------------------------
mdat= [];
[x,w]= gauleg(35,0,13);   % Gauss Legendre weights
x=flip(x);w=flip(w);      % proceed from small to large x
Nl=length(x);
la0 = 15;
z0=1.1
for i=1:Nl
  x(i)
  la=x(i)*la0;
  [z, fval, info] = fzero("f",z0);
  D = sum(na.*(aa*rha+ba*rha.^2));
  rha0=la/mu;
  D0 = sum(na)*(aa*rha0+ba*rha0^2);
  mdat = [mdat; x(i), z, D, D0, sum(na.*(rha==0)), sum(na.*(rha==1))];
  z0=z
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%		that's it
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function y = f (z)
global mu mua na rha la aa ba
  rha = min(max((z*mua-aa)/(2*ba),0),1);
  y = sum(na.*mua.*rha)-la;
end
%---------------------- Gauss Legendre Weight-------------------------
function [x,w]=gauleg(n, a, b)
if nargin == 1
    a = -1;
    b =  1;
end;
x = zeros(1, n);
w = zeros(1, n);
m = (n+1)/2;
h = b-a;

for ii=1:m
    z = cos(pi*(ii-.25)/(n+.5)); % Initial estimate.
    z1 = z+1;
    while abs(z-z1)>eps
        p1 = 1;
        p2 = 0;
        for jj = 1:n
            p3 = p2;
            p2 = p1;
            p1 = ((2*jj-1)*z*p2-(jj-1)*p3)/jj; % The Legendre polynomial.
        end
        pp = n*(z*p1-p2)/(z^2-1); % The L.P. derivative.
        z1 = z;
        z = z1-p1/pp;
        end
    x(ii) = z; % Build up the abscissas.
    x(n+1-ii) = -z;
    w(ii) = h/((1-z^2)*(pp^2)); % Build up the weights.
    w(n+1-ii) = w(ii);
        end

if a ~= -1 || b ~= 1
    x = (x+1)*(h/2) + a;
end
end
