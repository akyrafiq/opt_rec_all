%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%   otimize resource allocation for quadratic cost-function
%	z == Phi_0 Lagrange parameter controling sum of loads
%	Determine lambda as function of z =lammbda*phi0 !!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
global N M0 b c la z z0  mui mu
%--------------------------------------------------------------------
%  initialization
%--------------------------------------------------------------------
N  = 10000;              % number of processes
M0 =    1;		% exponentially distibituted resources MEAN
a  =    1;		% parameters of cost function
b  =    1;
c  =    1;
mui = rand(1,N);

SET = [1,2,5,10];
for p = SET;
mu0=1;
Beta=0.5*p;
mui = gamrnd(mu0,Beta,[1,N]);		% exprnd(M0,1,N);	% resource parameters
		% total capacity
%mui = exprnd(p,[1,N]);
mu  = sum(mui);
%--------------------------------------------------------------------
mdat= [];
zi= 0.0001:0.01:50;
Nl=length(zi);


for i=1:Nl
  z =zi(i);
  la = f(z);
  x_i = la/mu;
  D = N*(x_i.*(1+ 0.5*b*x_i +c*x_i.^2/3))/la;
  if ((la > 0) & (la < mu))
    phi0 = z/la;
    x0i = max(4*c*(z*mui-1),0);
    x0i = min((sqrt(b^2 +x0i) -b)/(2*c),1);			
    D_OPT = sum(x0i.*(1+ 0.5*b*x0i +c*x0i.^2/3))/la;
    
    mdat = [mdat; la/mu, phi0, sum((x0i <=0 )),sum((x0i >= 1)),D_OPT,D,D-D_OPT];
  end
end
scatter(mdat(:,1),mdat(:,3),'x')
hold on
plot(mdat(:,1),mdat(:,4),'o')
hold on
plot(mdat(:,1),mdat(:,5))
hold on
plot(mdat(:,1),mdat(:,6))
%plot(mdat(:,1),mdat(:,2))
end
set(gca, 'YScale', 'log');
xlim([0,1])
xlabel('Load ($\rho$)','Interpreter','latex')
%ylabel(' Langrange Parameter ($\phi_0$)','Interpreter','latex')
ylabel('Delay ($D$), Opmitmal Delay ($D^*$), No. of Idle ($N_0$) and Full ($N_F$) Resources ','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
function la = f (z)
global b c mu mui
  x0 = max(4*c*(mui*z-1),0);
  x0 = min((sqrt(b^2 +x0) -b)/(2*c),1);
  la = sum(mui.*x0);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%		that's it
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
