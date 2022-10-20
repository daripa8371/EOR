%% Solving Saturation Equations
% code to compute solution of saturation,concentration and 
% surfactant equations by Modified Method of Characteristics
% using explicit formulation (Yuan Yi-Rang 1993) and implicit finite
% difference method
%%
function [UU,CC,GG,ocut,wcut,ROIP] = nmmoc_surf_mod(u,v,S,C,G,miua,para,sigma,c0,g0,src)

format long

global miuo swr0 sor0 dt KK s0;

m=size(S,2); n=size(S,1); % Notice m = para.box.m+1 and n = para.box.n+1
dx = para.box.dx; dy = para.box.dy;
Q = S;
dtcal = dt*ones(n,m);   %make an array of dt
% define constant parameters for Pc
omega1 = 0.1; omega2 = 0.4;
phi = 1; %porosity
% fixed coordinates
[x,y] = meshgrid(para.box.left:para.box.dx:para.box.right,...
    para.box.bottom:para.box.dy:para.box.top);

% define critical capillary numbers 
% ie $$N_c $$ at which $$s_{ro}$$ and $$ s_{ra}$$ begin to decrease
Nco0  = 10^(-5);
Nca0 = 10^(-5);  %%% these two do not have to be the same

%% Parameter definitions
% Recompute residual saturations using (n+1)th time velocities and then 
% define normalized saturations of water and oil at IFT sigma as
% $$ \bar{s} = \frac{s-s_{ra}}{1-s_{ra}} $$
%
% $$ \tilde{s} = \frac{s-s_{ra}}{1-s_{ra}-s_{ro}} $$
[swr,sor] = compres(sigma,u,v,miua);
nsw = (Q-swr)/(1-swr);
nso = (Q-swr)/(1-swr-sor);

% recompute mobilities (with surfactants)
lambda_a = compmob(Q,C,sor,swr,1,1);
lambda_o = compmob(Q,C,sor,swr,0,1);
lambda = lambda_a + lambda_o;

% recompute fractional flow $$f = \lambda_a/\lambda $$  
% and $$ D = K(x) \lambda_o f $$
f = lambda_a./lambda;
D = KK.*lambda_o.*f;

% derivative of IFT wrt surf conc $$\sigma_\Gamma =-10.001/(\Gamma+1)^2$$
sigma_g = -10.001./(G+1).^2;

% compute capillary number
nca = sqrt(u^2+v^2).*miua./sigma; nco = sqrt(u^2+v^2).*miuo./sigma;
Nca = norm(nca); % compute 2-norm of Nc matrix ie largest singular value
Nco = norm(nco);

% recomputing derivative of residual saturations wrt surf conc
% $$ \frac{\partial s_{ra}}{\partial \Gamma}, \frac{\partial s_{ro}}{\partial \Gamma} $$
swr_g = zeros(n,m); sor_g =swr_g;

for j = 1:n  % columns
    for i = 1:m  % rows
        if Nca >= Nca0
            swr_g(j,i) = -(swr0*0.1534*10.001*Nca0^0.1534)/((sqrt(u(j,i)^2+v(j,i)^2)...
                *miua(j,i))^(0.1534)*sigma(j,i)^.8466*(G(j,i)+1)^2);
        end
        if Nco >= Nco0
            sor_g(j,i) = -(sor0*0.5213*10.001*Nco0^0.5213)/((sqrt(u(j,i)^2+v(j,i)^2)...
                *miua(j,i))^(0.5213)*sigma(j,i)^.4787*(G(j,i)+1)^2);
        end
    end
end

% derivatives of normalized saturations wrt surf conc
nsw_g = swr_g .*(Q-1)./(1-swr)^2;
nso_g = (swr_g.*(Q+sor-1)+sor_g.*(Q-swr))/(1-swr-sor)^2;

% derivative of relative perm wrt saturation
kra_s = 2.5*swr*(3*(nsw).^2-1)+1;
kro_s = 10*sor*nso - 5*sor-1;

% derivative of relative perm wrt surf conc
kra_g = 2.5*swr_g*(nsw.^3-nsw)+ (Q-1).*(2.5*swr*(3*nsw.^2-1)+ 1).*nsw_g/(1-swr)^2;
kro_g = 1-5*sor*nso + (1-nso).*(1-5*nso.*sor_g) - (1+5*sor-10*sor*nso).*nso_g;

% derivative of fractional flow func wrt saturation, pol conc and surf conc
f_s = kra_s.*lambda_o./(lambda.^2.*miua) - kro_s.*lambda_a./(lambda.^2*miuo);
%f_c = - (lambda_o.*lambda_a*miuo)./(lambda.^2.*miua); %to be redefined later
%f_g = kra_g.*lambda_o./(lambda.^2.*miua) - kro_g.*lambda_a./(lambda.^2*miuo);

% capillary pressure and its derivatives
pc = (sigma.*omega2*phi^.5)./(KK.^(.5)*(1-nso)^(1/omega1));
pc_s = pc./(omega1*(1-nso));
pc_g = (pc./sigma).*sigma_g + pc_s;


%% Solving for saturation
% Compute characteristics and use finite difference discretization for
% saturation equations on redefined coordinates.

% redefined coordinates
[xmod,ymod,dtmod] = eval_Xsurf(x,y,Q,Q,G,f,f_s,D,pc_s,pc_g,u,v,dtcal,para,1);

% bilinear interpolant for saturation on redefined coordinates
Qmod = interp2(x,y,Q,xmod,ymod);

% optimally redefined coordinates
[xmod,ymod,dtmod] = eval_X_surf(x,y,b_s,umod,vmod,dt);
% saturation values on redefined coordinates interpolated from knot values
Qmod = interp2(x,y,Q,xmod,ymod,'linear',0);


%redefined normalized saturations of water and oil at IFT sigma
%nsw = (Qmod-swr)/(1-swr);
nso = (Qmod-swr)/(1-swr-sor);

%updating coefficients with interpolated saturation( with surfactants )
lambda_a = compmob(Qmod,C,sor,swr,1,1);
lambda_o = compmob(Qmod,C,sor,swr,0,1);
lambda = lambda_a + lambda_o;
f = lambda_a./lambda;
D = KK.*lambda_o.*f;
pc_s = pc./(omega1*(1-nso));
pc_g = (pc./sigma).*sigma_g + pc_s;
f_c = - (lambda_o.*lambda_a*miuo)./(lambda.^2.*miua);
f_g = kra_g.*lambda_o./(lambda.^2.*miua) - kro_g.*lambda_a./(lambda.^2*miuo);

%intermediate calculation param
Ds = D.*pc_s;
Dg = D.* pc_g;

%------ test model for D ------
%D = 0.04*0.01.*Qmod.*(1-Qmod);
    
iter=1;
%preallocating for speed
AAA = zeros(n*m);
DDD = zeros(n*m,1);

while(iter<=(m)*(n-1)+1)
    cnt = (iter -1)/m; %% cnt = 0,1,2,3.... for iter = 1, m+1, 2m+1, 3m+1 etc...
    BB=zeros(n,m);AA=BB;CC=BB;DD=zeros(m,1);
    for i=1:1:m  %% i = x direction discretization formula and columns ie 2nd coordinate in matrix code
        for j=1:1:n %% j = y direction discretization formula and rows ie 1st coordinate in matrix code
            if j==i
                if iter == 1            % 1st or bottom row of grid
                    if i == 1           % 1st or left column
                        DD(i) = 1; %src
                        BB(j,i) = 1;
                    elseif i == m       % last or rightmost column of grid
                        DD(i) = Qmod(cnt+1,i);
                        BB(j,i) = 1;
                    else
                        DD(i) = Qmod(cnt+1,i);
                        BB(j,i) = 1;
                    end
                elseif iter == (m)*(n-1)+1   %topmost row of grid
                    if i == 1                % leftmost column
                        DD(i) = Qmod(cnt+1,i);
                        BB(j,i)= 1;
                    elseif i == m            % rightmost column
                        DD(i) = Qmod(cnt+1,i);
                        BB(j,i) = 1;
                    else
                        DD(i) = Qmod(cnt+1,i);
                        BB(j,i) = 1;
                    end
                else                         % interior rows of grid
                    if i == 1                % left most column of grid
                        DD(i) = Qmod(cnt+1,i);
                        BB(j,i) = 1;
                    elseif i == m            % right most column of grid
                        DD(i) = Qmod(cnt+1,i);
                        BB(j,i) = 1;
                    else
                        DD(i) =(1/dtmod(cnt+1,i))*Qmod(cnt+1,i) - f_c(cnt+1,i)*(u(cnt+1,i)*(C(cnt+1,i+1)-C(cnt+1,i-1))...
                            /(2*dx)+v(cnt+1,i)*(C(cnt+2,i)-C(cnt,i))/(2*dy))-f_g(cnt+1,i)*(u(cnt+1,i)*...
                            (G(cnt+1,i+1)-G(cnt+1,i-1))/(2*dx)+v(cnt+1,i)*(G(cnt+2,i)-G(cnt,i))/(2*dy))...
                            -(Dg(cnt+1,i+1)/(2*dx*dx)*(G(cnt+1,i+1)-G(cnt+1,i)) + Dg(cnt+1,i-1)/(2*dx*dx)...
                            *(G(cnt+1,i-1)-G(cnt+1,i)) + Dg(cnt+1,i)/(2*dx*dx)*(G(cnt+1,i-1)+G(cnt+1,i+1)...
                            -2*G(cnt+1,i)) + Dg(cnt+2,i)/(2*dy*dy)*(G(cnt+2,i)-G(cnt+1,i)) + Dg(cnt,i)/...
                            (2*dy*dy)*(G(cnt,i)-G(cnt+1,i)) + Dg(cnt+1,i)/(2*dy*dy)*(G(cnt+2,i)+G(cnt,i)...
                            -2*G(cnt+1,i)));
                        AA(j,i)   =  (Ds(cnt,i)+Ds(cnt+1,i))/(2*dy*dy);
                        CC(j,i)   =  (Ds(cnt+1,i)+Ds(cnt+2,i))/(2*dy*dy);
                        BB(j,i)   = 1/dtmod(cnt+1,i)-((1/(2*dx*dx))*(Ds(cnt+1,i-1)+2*Ds(cnt+1,i)+Ds(cnt+1,i+1))+...
                            (1/(2*dy*dy))*(Ds(cnt,i)+2*Ds(cnt+1,i)+Ds(cnt+2,i)));
                        BB(j,i+1) =  (Ds(cnt+1,i)+Ds(cnt+1,i+1))/(2*dx*dx);
                        BB(j,i-1) =  (Ds(cnt+1,i-1)+Ds(cnt+1,i))/(2*dx*dx);
                    end
                end      
%                 if iter ~= 1 && iter ~= m*(n-1)+1                         % interior rows of grid
%                     if i ~= 1 && i ~= m
%                         DD(i) =(1/dtmod(cnt+1,i))*Qmod(cnt+1,i) - f_c(cnt+1,i)*(u(cnt+1,i)*(C(cnt+1,i+1)-C(cnt+1,i-1))...
%                             /(2*dx)+v(cnt+1,i)*(C(cnt+2,i)-C(cnt,i))/(2*dy))-f_g(cnt+1,i)*(u(cnt+1,i)*...
%                             (G(cnt+1,i+1)-G(cnt+1,i-1))/(2*dx)+v(cnt+1,i)*(G(cnt+2,i)-G(cnt,i))/(2*dy))...
%                             -(Dg(cnt+1,i+1)/(2*dx*dx)*(G(cnt+1,i+1)-G(cnt+1,i)) + Dg(cnt+1,i-1)/(2*dx*dx)...
%                             *(G(cnt+1,i-1)-G(cnt+1,i)) + Dg(cnt+1,i)/(2*dx*dx)*(G(cnt+1,i-1)+G(cnt+1,i+1)...
%                             -2*G(cnt+1,i)) + Dg(cnt+2,i)/(2*dy*dy)*(G(cnt+2,i)-G(cnt+1,i)) + Dg(cnt,i)/...
%                             (2*dy*dy)*(G(cnt,i)-G(cnt+1,i)) + Dg(cnt+1,i)/(2*dy*dy)*(G(cnt+2,i)+G(cnt,i)...
%                             -2*G(cnt+1,i)));
%                         AA(j,i)   =  (Ds(cnt+1-1,i)+Ds(cnt+1,i))/(2*dy*dy);
%                         CC(j,i)   =  (Ds(cnt+1,i)+Ds(cnt+1+1,i))/(2*dy*dy);
%                         BB(j,i)   = 1/dtmod(cnt+1,i)-((1/(2*dx*dx))*(Ds(cnt+1,i-1)+2*Ds(cnt+1,i)+Ds(cnt+1,i+1))+...
%                             (1/(2*dy*dy))*(Ds(cnt+1-1,i)+2*Ds(cnt+1,i)+Ds(cnt+1+1,i)));
%                         BB(j,i+1) =  (Ds(cnt+1,i)+Ds(cnt+1,i+1))/(2*dx*dx);
%                         BB(j,i-1) =  (Ds(cnt+1,i-1)+Ds(cnt+1,i))/(2*dx*dx);
%                     end
%                 end              
                
            end
        end
    end
    if cnt == 0
        AAA(1:n,1:2*m) = [BB CC];
    elseif cnt == n-1
        AAA((m-1)*n+1:m*n,(n-2)*m+1:n*m) = [AA BB];        
    else
        AAA(cnt*n+1:(cnt+1)*n,(cnt-1)*m+1:(cnt+2)*m) = [AA BB CC];
    end
    DDD(cnt*m+1:(cnt+1)*m) = DD;   
    iter = iter+m;
end

% QQQ = AAA\DDD;
% Qnew = reshape(QQQ,m,n)';

Qnew=(reshape(bicgstab(AAA,DDD,1*10^(-10),600,[],[]),m,n))';

% imposing homogeneous neumann boundary conditions 

for i=1:1:m  %% i = x direction discretization formula and columns ie 2nd coordinate in matrix code
        for j=1:1:n %% j = y direction discretization formula and rows ie 1st coordinate in matrix code
            
            if j == 1            % 1st or bottom row of grid
                if i == 1           % 1st or left column
                    Qnew(j,i) =1;
                elseif i == m       % last or rightmost column of grid
                    Qnew(j,i) = Qnew(j,i-1)/2+Qnew(j+1,i)/2; %Qnew(j,i-1);
                end
            elseif j == n   %topmost row of grid
                if i == 1                % leftmost column
                    Qnew(j,i) = Qnew(j,i+1)/2 + Qnew(j-1,i)/2; %Qnew(j-1,i)/2; % + Qnew(j,i)/2; %
                elseif i == m            % rightmost column
                    Qnew(j,i) = Qnew(j,i-1)/2 + Qnew(j-1,i)/2;
                end
            else                         % interior rows of grid
                if i == 1                % left most column of grid
                    Qnew(j,i) = (2*Qnew(j,i+1)+Qnew(j+1,i)+Qnew(j-1,i))/4; %(Qnew(j,i+1)+Qnew(j,i)+Qnew(j-1,i))/3; %
                elseif i == m            % right most column of grid
                    Qnew(j,i) = (2*Qnew(j,i-1)+Qnew(j+1,i)+Qnew(j-1,i))/4; %(Qnew(j,i-1)+Qnew(j,i)+Qnew(j-1,i))/3; %
                end
            end
            
        end
end


Qnew(Qnew>1) = 1;
Smax=max(max(Qnew)); disp(Smax);
%% Solving for concentration of polymer
% recompute characteristics for concentration equation and discretize using
% finite difference on recomputed coordinates

%new redefined coordinates for concentration equation
[xmod2,ymod2,dtmod] = eval_Xsurf(x,y,Q,Qnew,G,f,f_s,D,pc_s,pc_g,u,v,dtcal,para,2);

%bilinear interpolant for saturation on redefined coordinates
Cmod = interp2(x,y,C,xmod2,ymod2);
%Cnew=zeros(n,m);  % initialize Cnew for speed of computation
iter =1;
AAA = zeros(n*m);
DDD = zeros(n*m,1);

while(iter<=(m)*(n-1)+1)
    cnt = (iter -1)/m; %% cnt = 0,1,2,3.... for iter = 1, m+1, 2m+1, 3m+1 etc...
    BB=zeros(n,m);AA=BB;CC=BB;DD=zeros(m,1);
    for i=1:1:m  %% i = x direction discretization formula and columns ie 2nd coordinate in matrix code
        for j=1:1:n %% j = y direction discretization formula and rows ie 1st coordinate in matrix code
            if j==i
                if iter == 1
                    if i == 1
                        DD(i) = c0; %*src/Qnew(cnt+1,i);
                        BB(j,i) =1;%1/dtmod(cnt+1,i)+1/Qnew(cnt+1,i);

                    elseif i == m
                        DD(i) = Cmod(cnt+1,i);
                        BB(j,i) = 1;

                    else
                        DD(i) = Cmod(cnt+1,i);
                        BB(j,i) = 1;
                    end
                elseif iter == (m)*(n-1)+1
                    if i == 1
                        DD(i) = Cmod(cnt+1,i);
                        BB(j,i)= 1;

                    elseif i == m
                        DD(i) = Cmod(cnt+1,i);
                        BB(j,i) = 1;
                    else
                        DD(i) = Cmod(cnt+1,i);
                        BB(j,i) = 1;
                    end
                else
                    if i == 1
                        DD(i) = Cmod(cnt+1,i);
                        BB(j,i) = 1;
                    elseif i == m
                        DD(i) = Cmod(cnt+1,i);
                        BB(j,i) = 1;
                    else                
                        DD(i) = (1/dtmod(cnt+1,i))*Cmod(cnt+1,i);
                        BB(j,i) =  1/dtmod(cnt+1,i);
                    end
                end             
            end
        end
    end
    if cnt == 0
        AAA(1:n,1:2*m) = [BB CC];
    elseif cnt == n-1
        AAA((m-1)*n+1:m*n,(n-2)*m+1:n*m) = [AA BB];        
    else
        AAA(cnt*n+1:(cnt+1)*n,(cnt-1)*m+1:(cnt+2)*m) = [AA BB CC];
    end
    DDD(cnt*m+1:(cnt+1)*m) = DD; 
    iter = iter+m;
end

Cnew=(reshape(bicgstab(AAA,DDD,1*10^(-10),600,[],[]),m,n))';

% imposing homogeneous neumann boundary conditions 

for i=1:1:m  %% i = x direction discretization formula and columns ie 2nd coordinate in matrix code
        for j=1:1:n %% j = y direction discretization formula and rows ie 1st coordinate in matrix code
            
            if j == 1            % 1st or bottom row of grid
                if i == 1           % 1st or left column
                    Cnew(j,i) =c0;
                elseif i == m       % last or rightmost column of grid
                    Cnew(j,i) = Cnew(j,i-1)/2+Cnew(j+1,i)/2;
                end
            elseif j == n   %topmost row of grid
                if i == 1                % leftmost column
                    Cnew(j,i) = Cnew(j,i+1)/2 + Cnew(j-1,i)/2;
                elseif i == m            % rightmost column
                    Cnew(j,i) = Cnew(j,i-1)/2 + Cnew(j-1,i)/2;
                end
            else                         % interior rows of grid
                if i == 1                % left most column of grid
                    Cnew(j,i) = (2*Cnew(j,i+1)+Cnew(j+1,i)+Cnew(j-1,i))/4;
                elseif i == m            % right most column of grid
                    Cnew(j,i) = (2*Cnew(j,i-1)+Cnew(j+1,i)+Cnew(j-1,i))/4;
                end
            end
            
        end
end

Cnew(Cnew>c0) = c0;
Cmax=max(max(Cnew)); disp(Cmax);
%% Solving for Surfactant concentration
% recompute characteristics for surfactant equation and discretize using
% finite difference on recomputed coordinates
     
%new redefined coordinates for concentration equation
[xmod2,ymod2,dtmod] =  eval_Xsurf(x,y,Q,Qnew,G,f,f_s,D,pc_s,pc_g,u,v,dtcal,para,3);

%bilinear interpolant for sur conc on redefined coordinates
Gmod = interp2(x,y,G,xmod2,ymod2);

%updating coefficients with interpolated Surf conc 
sigmamod = 10.001./(Gmod+1)-0.001;
sigma_g_mod = -10.001./(Gmod+1).^2;
[swr,sor] = compres(sigmamod,u,v,miua);
nso = (Qmod-swr)/(1-swr-sor);
lambda_a = compmob(Qmod,C,sor,swr,1,1);
lambda_o = compmob(Qmod,C,sor,swr,0,1);
lambda = lambda_a + lambda_o;
f = lambda_a./lambda;
D = KK.*lambda_o.*f;
pc_s = pc./(omega1*(1-nso));
pc_g = (pc./sigmamod).*sigma_g_mod + pc_s;

%Intermediate parameter for code
F = D.*pc_g./Qnew;
iter =1;
AAA = zeros(n*m);
DDD = zeros(n*m,1); 

while(iter<=(m)*(n-1)+1)
    cnt = (iter -1)/m; %% cnt = 0,1,2,3.... for iter = 1, m+1, 2m+1, 3m+1 etc...
    BB=zeros(n,m);AA=BB;CC=BB;DD=zeros(m,1);
    for i=1:1:m  %% i = x direction discretization formula and columns ie 2nd coordinate in matrix code
        for j=1:1:n %% j = y direction discretization formula and rows ie 1st coordinate in matrix code
            if j==i
                if iter == 1
                    if i == 1
                        DD(i) = g0;%/Qnew(cnt+1,i);
                        BB(j,i) =1;%1/dtmod(cnt+1,i) + 1/Qnew(cnt+1,i);
                    elseif i == m
                        DD(i) = Gmod(cnt+1,i);
                        BB(j,i) = 1;
                    else
                        DD(i) = Gmod(cnt+1,i);
                        BB(j,i) = 1;
                    end
                elseif iter == (m)*(n-1)+1
                    if i == 1
                        DD(i) = Gmod(cnt+1,i);
                        BB(j,i)= 1;
                    elseif i == m
                        DD(i) = Gmod(cnt+1,i);
                        BB(j,i) = 1;
                    else
                        DD(i) = Gmod(cnt+1,i);
                        BB(j,i) = 1;
                    end
                else
                    if i == 1
                        DD(i) = Gmod(cnt+1,i);
                        BB(j,i) = 1;
                    elseif i == m
                        DD(i) = Gmod(cnt+1,i);
                        BB(j,i) = 1;
                    else
                        
                        DD(i) = (1/dtmod(cnt+1,i))*Gmod(cnt+1,i);
                        AA(j,i)   =  F(cnt+1,i)/(dy*dy);
                        CC(j,i)   =  F(cnt+1,i)/(dy*dy);
                        BB(j,i)   =  1/dtmod(cnt+1,i) - 2*(1/dx^2+1/dy^2)*F(cnt+1,i);
                        BB(j,i+1) =  F(cnt+1,i)/(dx*dx);
                        BB(j,i-1) =  F(cnt+1,i)/(dx*dx);
                    end
                end
                
            end
        end
    end
    if cnt == 0
        AAA(1:n,1:2*m) = [BB CC];
    elseif cnt == n-1
        AAA((m-1)*n+1:m*n,(n-2)*m+1:n*m) = [AA BB];        
    else
        AAA(cnt*n+1:(cnt+1)*n,(cnt-1)*m+1:(cnt+2)*m) = [AA BB CC];
    end
    DDD(cnt*m+1:(cnt+1)*m) = DD; 

    iter = iter+m;
end

Gnew=(reshape(bicgstab(AAA,DDD,1*10^(-10),600,[],[]),m,n))';
% imposing homogeneous neumann boundary conditions 

for i=1:1:m  %% i = x direction discretization formula and columns ie 2nd coordinate in matrix code
        for j=1:1:n %% j = y direction discretization formula and rows ie 1st coordinate in matrix code
            
            if j == 1            % 1st or bottom row of grid
                if i == 1           % 1st or left column
                    Gnew(j,i) =g0;
                elseif i == m       % last or rightmost column of grid
                    Gnew(j,i) = Gnew(j,i-1)/2+Gnew(j+1,i)/2;
                end
            elseif j == n   %topmost row of grid
                if i == 1                % leftmost column
                    Gnew(j,i) = Gnew(j,i+1)/2 + Gnew(j-1,i)/2;
                elseif i == m            % rightmost column
                    Gnew(j,i) = Gnew(j,i-1)/2 + Gnew(j-1,i)/2;
                end
            else                         % interior rows of grid
                if i == 1                % left most column of grid
                    Gnew(j,i) = (2*Gnew(j,i+1)+Gnew(j+1,i)+Gnew(j-1,i))/4;
                elseif i == m            % right most column of grid
                    Gnew(j,i) = (2*Gnew(j,i-1)+Gnew(j+1,i)+Gnew(j-1,i))/4;
                end
            end
            
        end
end
if g0 ~= 0
Gnew(Gnew>g0) = g0; Gmax =max((max(Gnew))); disp(Gmax);
end
disp('mmoc calculation done');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Visualization studies

% % figure(100)
% % plot(x(1,:),S(2,:),'Color','blue');
% % hold on
% % axis([0 1 -0.5 0.9])
% % plot(x(1,:),Qnew(2,:)-0.4,'Color','red');
% 
% 
% % figure(101)
% % plot(x(1,:),S(floor(n-2),:),'Color','blue');
% % hold on
% % axis([0 1 -0.5 0.9])
% % plot(x(1,:),Qnew(floor(n-2),:)-0.4,'Color','red');


% pause

UU=Qnew;
CC=Cnew;
GG=Gnew;
%%%% Oil Recovered calculations 
ocut = lambda_o(n,m)*src./lambda(n,m);
wcut = lambda_a(n,m)*src*Qnew(n,m)./lambda(n,m);
ROIP = sum(sum(1 - Qnew));



