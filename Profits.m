%% Generating Servers%%
N=100;
mui = rand(1,N);
for p=1:N
mui2(p)=mui(p)^2;
end
mu=sum(mui);
%%%%%%%%%%%%%%%%%%%%%%%
for i=1:40
%% Distribution Parameters %%
l_bar=20;
l0=N*l_bar;
alpha=20;
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Revenue Parameters %%
l_max=i;
P=2;
r=10;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cost Function Parameters %%
zeta=1;
Zeta0=.2;
Zeta1=.2;
Zeta2=.2;
A=5;
A0=1;
A1=1;
A2=1;
for k=1:N
    exmui(k)= exp(A1*mui(k));
end
exmu=sum(exmui);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Loss Function %%
L=1;
L0=1;
L1=1;
L2=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
%% Obtaining Phi0%%
 for k = 1:1000
lambda = mu*k/1000;
beta0 = fzero(@(beta) Sfun(beta,mui,mu,lambda), [0 200]);
phi0 = 1/(lambda*beta0);
chek = sum( sqrt(mui/(lambda*phi0)).*heaviside(mui-sqrt(mui/(lambda*phi0))) ... 
        + mui.*heaviside(sqrt(mui/(lambda*phi0))-mui)) - (mu - lambda);
max(abs(chek));  % should be small
Phi(k)=phi0;
LOAD(k)= lambda;
 end
 %% Optimal Profits %%
 Poo= sum(Wfunc(LOAD(k)/l0,alpha)*(sqrt(phi0.*mui)-1).*heaviside(sqrt(phi0.*mui)-1));
 P_Opt(k)=r*l0*alpha -sum(Poo);
 %}
 %% Full Penalty%%
P_FULL= P*(l0*alpha - l_max);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Full Revenue%%10
R = l0*(r+P)*gammainc(alpha+1,l_max)-l_max*gammainc(alpha,l_max);
%%%%%%%%%%%%%%%%
%% Expected Profits %%
Lin_PFT(i) = (R-(N*zeta*l0)/mu)*gammainc(alpha+1,l_max)-P_FULL-A*sum(mui);
Quad_PFT(i) = (R-l_max*gammainc(alpha,l_max) - N*l0*(Zeta0*gammainc(alpha,l_max)/l0...
    + Zeta1*gammainc(alpha+1,l_max)/mu + Zeta2*l0*gammainc(alpha+2,l_max)/(mu^2)))- P_FULL-(N*A0+A1*mu+A2*sum(mui2));
%Exp_PFT(i) = (R-(A0*exmu-l_max)*gammainc(alpha,l_max))-N*Zeta0(1/(1+(Zeta1*l0)/mu))^alpha*gammainc(alpha,l_max)-P_FULL;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Optimal Profits %%
%P(i) = r*l0*alpha - sum(integral(sqrt(mui*phi0)-1).*heaviside(sqrt(mui*phi0)-1),-Inf,Inf);
%%%%%%%%%%%%%%%%%%%
X_VAL(i)= l_max;
end
%% Graphic Reuslts%%
plot(X_VAL,Lin_PFT)
hold on
plot (X_VAL,Quad_PFT)
%set(gca, 'YScale', 'log');
xlabel('Maximum Load','Interpreter','latex') 
ylabel('Profits','Interpreter','latex') 
set(gca,'TickLabelInterpreter','latex')
%% Function to obtain phi0 %%
function S = Sfun(beta,mui,mu,lambda)
  S = sum( sqrt(beta*mui).*heaviside(mui-beta) ...
      + mui.*heaviside(beta-mui)) - (mu - lambda);
end
%% Weight Function
function Weight=Wfunc(x,y);
Weight = x^y*exp(-x)/gamma(y+1);
end