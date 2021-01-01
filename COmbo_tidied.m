 
    %% Start %%
    clear;clc;
    mdat_ad=[];
    mdat_opt=[];
    l_0=15;
    alpha=2
    for v=4:0.25:10
    %% Setting Maximum Load and Obtaining the servers%%
    x_max=v;
    U=5;
    U1=U^(0);
    U2=U;
    U3=U^2; %budget paramters-must be sublinear
    Powers=[U1,U2,U3];
    Conf=[];
    lambda_max=x_max*l_0;
    for i=0:floor(lambda_max/U1)
        for j= 0:floor(lambda_max/U2)
            for k=0:floor(lambda_max/U3)
                    if i*U1+U2*j+U3*k>=lambda_max & i*U1+U2*j+U3*k<=1.2*lambda_max; %Checking if combination is feasible
                        Config = [i,j,k];
                        N= i*U1+U2*j+U3*k; %Power of the configuration
                        Conf=[Conf;Config N];
                    else ;
                    end
                end
            end
    end
    Combos = size(Conf,1);
    %% Configuring the budget and Baseline Costs %%
    B=4; %Budget Parameters
    B1=B^0;
    B2=B^1;
    B3=B^2;
    A=7; %Baseline Parameters
    A1=A^0;
    A2=A;
    A3=A^2;
    Costs=[];
    for i=1:length(Conf);
        B_mu= B1*Conf(i,1) + B2*Conf(i,2)+B3*Conf(i,3);
        C_mu= A1*Conf(i,1) + A2*Conf(i,2)+A3*Conf(i,3);
        Costs=[Costs;B_mu, C_mu];
    end
    Conf_wB=[Conf,Costs];
    %% Running Costs %%
    aa=1;
    c_a=2;
    bb=0.5;
    c_1=[1,2,4];
    C1=[];
    for y=1:length(Conf_wB);
    mui = [U1*ones(1,Conf_wB(y,1)),U2*ones(1,Conf_wB(y,2)),U3*ones(1,Conf_wB(y,3))];
    Costs = [c_1(1)*ones(1,Conf_wB(y,1)),c_1(2)*ones(1,Conf_wB(y,2)),c_1(3)*ones(1,Conf_wB(y,3))];
    mu=Conf_wB(y,4);
    rho_ad= [U1*lambda_max/mu, U2*lambda_max/mu, U3*lambda_max/mu];
    C1_ad=[];
    for t=1:3;
        %C1_ad(t)=c_1(t)*(rho_ad(t)+0.5*rho_ad(t)^2)*Conf_wB(y,t);
        C1_ad(t)=l_0*c_1(t)*(Powers(t)*aa*gammainc(alpha+1,x_max)/mu + Powers(t)^2*bb*l_0*gammainc(alpha+2,x_max)/(mu^2))*Conf_wB(y,t);
    end
    C1_adt=sum(C1_ad);
    % This seciton here regaridng optimal reosurces needs editing
    GL_Weights= gauleg(length(mui), 0, x_max); %Reimer only sets N=
    for h=1:length(mui)
        z=x_max;
       L_Op(h) = max(4*Costs(h)*(z*mui(h)-1),0);
        L_Op(h) = min((sqrt(0.5^2 +L_Op(h)) -0.5)/(2*Costs(h)),1);
    end
    for t=1:3;
        C1_ad(t)=c_1(t)*(rho_ad(t)+0.5*rho_ad(t)^2)*Conf_wB(y,t);
    end
    C_opt=sum(L_Op);
    %End of Edited Section
    C1=[C1;C1_adt,C_opt];
    end
    Conf_wR=[Conf_wB, C1];
    %% Find Our optimal Solution %%
    C_Min_Ad= Conf_wR(1,5)+Conf_wR(1,6)+Conf_wR(1,7);
    C_Min_Opt= Conf_wR(1,5)+Conf_wR(1,6)+Conf_wR(1,8);
    OptComb_ad=1;
    OptComb_opt=1;
    for y=2:length(Conf_wR);
        if Conf_wR(y,5)+Conf_wR(y,6)+Conf_wR(y,7) < C_Min_Ad;
            OptComb_ad=y;
            C_Min_Ad = Conf_wR(y,5)+Conf_wR(y,6)+Conf_wR(y,7);
        end
    end
    for y=2:length(Conf_wR);
        if Conf_wR(y,5)+Conf_wR(y,6)+Conf_wR(y,8) < C_Min_Opt;
            OptComb_opt=y;
            C_Min_Opt = Conf_wR(y,5)+Conf_wR(y,6)+Conf_wR(y,8);
        end
    end
    Solution_ad=[Conf_wR(OptComb_ad,:), C_Min_Ad];
    Solution_opt=[Conf_wR(OptComb_opt,:), C_Min_Ad];
    %% Final Table %%
    mdat_ad=[mdat_ad;x_max, Solution_ad];
    mdat_opt=[mdat_opt;x_max, Solution_opt];
    end
    %% Tidying Results
    mdat_ad
    mdat_opt
    %lab=['X_Max';'N_1';'N_5';'N_25';'mu';'B-mu';'C-mu';'C1_ah';'Total Cost'];
    %Results=[lab;mdat]
     %% Visualise Reults
    plot(mdat_opt(:,1),mdat_opt(:,10))
    hold on
    plot(mdat_ad(:,1),mdat_ad(:,10))
    %% Optimal Allocation Solution
    function la = f (z)
    global b c mu mui
      x0 = max(4*c*(mui*z-1),0);
      x0 = min((sqrt(b^2 +x0) -b)/(2*c),1);
    end
 
    %% Gaus Leg
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
    %% P-Distribution of X
    function y= P_x(x)
    y=x*exp(-x)
    end
    %% Notes
     %{
    zi= 1:0.1:50;
    Nl=length(zi);
    for i=1:Nl;
      z =zi(i);
      la = f(z);
      if ((la > 0) & (lambda_max < mu));
        phi0 = z/lambda_max;
        x0i = max(4*c*(z*mui-1),0);
        x0i = min((sqrt(b^2 +x0i) -b)/(2*c),1);        
        C1 = sum(x0i.*(1+ 0.5*b*x0i +c*x0i.^2/3))/lambda_max;
      end
    end
    %}
   %% 
    function y=ff(z)
    rha = min(max((z*mua./c1a-aa)/(2*ba),0),1);
  y = sum(na.*rha.*mua)-la;
    end