%% Solving Saturation Equations
% code to compute solution of saturation,concentration and 
% surfactant equations by Modified Method of Characteristics
% using explicit formulation (Yuan Yi-Rang 1993) and implicit finite
% difference method
% This code implements Corey permeability without surfactant and Amaefule
% type with surfactant
%%
function [UU,CC,GG,ocut,wcut,ROIP] = nmmoc_surf_mod_neumann(u,v,S,C,G,miua,para,sigma,src)
% new formulation with 2nd order Neumann boundary conditions

format long

global miuo swr0 sor0 dt KK c0 g0 theta;
G1 = src; %source for water saturation equation
G2 = c0*src; %source for concentration equation
G3 = g0*src; %source for concentration equation
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
Nco0 = 1.44*10^(-4);  %% Values from Amafuele Handy 1982
Nca0 = 1.44*10^(-4);  %%% these two do not have to be the same

%% Parameter definitions
% Recompute residual saturations using (n+1)th time velocities and then 
% define normalized saturations of water and oil at IFT sigma as
% $$ \bar{s} = \frac{s-s_{ra}}{1-s_{ra}} $$
%
% $$ \tilde{s} = \frac{s-s_{ra}}{1-s_{ra}-s_{ro}} $$
[swr,sor] = compres(sigma,u,v,miua);


% recompute mobilities (with surfactants)
lambda_a = compmob(Q,C,sor,swr,1);
lambda_o = compmob(Q,C,sor,swr,0);
lambda = lambda_a + lambda_o;

% recompute fractional flow $$f = \lambda_a/\lambda $$  
% and $$ D = K(x) \lambda_o f $$
f = lambda_a./lambda;
D = KK.*lambda_o.*f;

% recomputing derivative of residual saturations wrt surf conc
% $$ \frac{\partial s_{ra}}{\partial \Gamma}, \frac{\partial s_{ro}}{\partial \Gamma} $$
swr_g = zeros(n,m); sor_g =swr_g; sigma_g = swr_g;

if g0 ~= 0
    % derivative of IFT wrt surf conc $$\sigma_\Gamma =-10.001/(\Gamma+1)^2$$
    sigma_g = -10.001./(G+1).^2; % == 0 when no surfactant
    nsw = (Q-swr)/(1-swr);
    nso = (Q-swr)/(1-swr-sor);
%     disp(nso(1,1))
    % compute capillary number
    nca = sqrt(u.^2+v.^2).*miua./sigma; nco = sqrt(u.^2+v.^2).*miuo./sigma;
    Nca = mean(mean(nca)); % compute 2-norm of Nc matrix ie largest singular value
    Nco = mean(mean(nco));
    nca_interm = mean(mean(sqrt(u.^2+v.^2))); nco_interm = mean(mean(sqrt(u.^2+v.^2).*miuo));
    for j = 1:n  % columns
        for i = 1:m  % rows
            if Nca >= Nca0
%                 swr_g(j,i) = -(swr0*0.1534*10.001*Nca0^0.1534)/((sqrt(u(j,i)^2+v(j,i)^2)...
%                     *miua(j,i))^(0.1534)*sigma(j,i)^.8466*(G(j,i)+1)^2);
                swr_g(j,i) = -(swr0*0.1534*10.001)*(Nca0/nca_interm/miua(j,i))^(0.1534)  ...
                    /sigma(j,i)^(.8466)/( G(j,i)+1 )^2 ;
            end
            if Nco >= Nco0
                sor_g(j,i) = -(sor0*0.5213*10.001)*(Nco0/nco_interm)^(0.5213)  ...
                    /sigma(j,i)^(.4787)/( G(j,i)+1 )^2 ;
%                 sor_g(j,i) = -(sor0*0.5213*10.001*Nco0^0.5213)/((sqrt(u(j,i)^2+v(j,i)^2)...
%                     *miua(j,i))^(0.5213)*sigma(j,i)^.4787*(G(j,i)+1)^2);
            end
        end
    end
    % derivative of relative perm wrt saturation
    kra_s = 2.5*swr*(3*(nsw).^2-1)+1;
    kro_s = 10*sor*nso - 5*sor-1;
    
    % derivatives of normalized saturations wrt surf conc
    nsw_g = swr_g .*(Q-1)./(1-swr)^2;  % == 0 when no surfactant
    nso_g = (swr_g.*(Q+sor-1)+sor_g.*(Q-swr))/(1-swr-sor)^2; % == 0 when no surfactant
    
    %
    % derivative of relative perm wrt surf conc
    kra_g = 2.5*swr_g*(nsw.^3-nsw)+ (Q-1).*(2.5*swr*(3*nsw.^2-1)+ 1).*nsw_g/(1-swr)^2; % == 0 when no surfactant
    kro_g =  -(1-nso).*(5*nso.*sor_g) - (1+5*sor-10*sor*nso).*nso_g; % == 0 when no surfactant
    
    pc = (omega2.*sigma.*sqrt(phi./KK))./((1-nso)^(1/omega1));
    pc_s = pc./(omega1*(1-nso));
    pc_g = (pc./sigma).*sigma_g + pc_s.*nso_g; % == 0 when no surfactant
else
    
    nsw = (Q-swr)/(1-swr-sor);
    %nsw^(-1/omega1)
    %pause
    kra_s = ((2+3*theta)/theta)*nsw.^((2+2*theta)/theta);
    kro_s = -2*(1-nsw).*(1-nsw.^((2+theta)/theta)) - ((2+theta)/theta)*(1-nsw).^2.*nsw^(2/theta);
    %kro_s = -2*(1-nso).*(1-nso.^(1.5))-1.5*((1-nso).^2).*sqrt(nso);
    
    % derivatives of normalized saturations wrt surf conc
    nsw_g = (swr_g .*(Q-1+sor)+sor_g.*(Q-swr))./(1-swr-sor)^2;  % == 0 when no surfactant
    
    kra_g = kra_s.*nsw_g;
    kro_g = kro_s.*nsw_g;
    
    % capillary pressure and its derivatives
    %pc = (sigma*omega2*sqrt(phi))./(sqrt(KK)*(1-nsw)^(1/omega1)); % Older form
    %pc_s = pc./(omega1*(1-nsw)); % Older form
    pc = ((omega2*sigma.*sqrt(phi./KK))./(nsw)^(1/omega1));
    pc_s = - pc./(omega1*nsw);
    pc_g = (pc./sigma).*sigma_g + pc_s.*nsw_g; % == 0 when no surfactant
end




% derivative of fractional flow func wrt saturation, conc and surf conc
f_s = kra_s.*lambda_o./(lambda.^2.*miua) - kro_s.*lambda_a./(lambda.^2*miuo);
%f_c = - (lambda_o.*lambda_a*miuo)./(lambda.^2.*miua); %to be redefined later
%f_g = kra_g.*lambda_o./(lambda.^2.*miua) - kro_g.*lambda_a./(lambda.^2*miuo);




%% Solving for saturation
% Compute characteristics and use finite difference discretization for
% saturation equations on redefined coordinates.

% redefined coordinates
[xmod,ymod] = eval_Xsurf_neumann(x,y,Q,Q,G,f,f_s,D,pc_s,pc_g,u,v,dtcal,para,1);


% fprintf('xmod(1,1) = %12.10f \nymod(1,1) = %12.10f \n',xmod(1,1),ymod(1,1))
% fprintf('xmod(n,m) = %12.10f \nymod(n,m) = %12.10f \n',xmod(n,m),ymod(n,m))



% bilinear interpolant for saturation on redefined coordinates
Qmod = interp2(x,y,Q,xmod,ymod);


% % optimally redefined coordinates
%  [xmod,ymod] = eval_X(x,y,b_s,umod,vmod,dt);
% % saturation values on redefined coordinates interpolated from knot values
% Qmod = interp2(x,y,Q,xmod,ymod,'linear',0);


%redefined normalized saturations of water and oil at IFT sigma
if g0~=0
    nsw = (Qmod-swr)/(1-swr); nso = (Qmod-swr)/(1-swr-sor);
    nso_g = (swr_g.*(Qmod+sor-1)+sor_g.*(Qmod-swr))/(1-swr-sor)^2;
else
    nsw = (Qmod-swr)/(1-swr-sor);
end
%nso = (Qmod-swr)/(1-swr-sor);

%updating coefficients with interpolated saturation( with surfactants )
lambda_a = compmob(Qmod,C,sor,swr,1);
lambda_o = compmob(Qmod,C,sor,swr,0);
lambda = lambda_a + lambda_o;
f = lambda_a./lambda;
f_c = - (lambda_o.*lambda_a*miuo)./(lambda.^2.*miua);
f_g = kra_g.*lambda_o./(lambda.^2.*miua) - kro_g.*lambda_a./(lambda.^2*miuo); % == 0 when no surfactant
D = KK.*lambda_o.*f;  % $$computing \bar{D} = D(\bar{s}^n, c^n)$$

if g0~=0
    pc = (omega2.*sigma.*sqrt(phi./KK))./((1-nso)^(1/omega1));
    pc_s = pc./(omega1*(1-nso));
    pc_g = (pc./sigma).*sigma_g + pc_s.*nso_g; % == 0 when no surfactant  
else
    pc = ((sigma*omega2*sqrt(phi))./(sqrt(KK))./(nsw)^(1/omega1));
    pc_s = - pc./(omega1*nsw); % $$computing \frac{\partial p_c}{\partial s}$$
    pc_g = (pc./sigma).*sigma_g + pc_s.*nsw_g; % $$computing \frac{\partial p_c}{\partial \Gamma}$$
end


%intermediate calculation param
Ds = D.*pc_s;
Dg = D.*pc_g; % == 0 when no surfactant

%------ test model for D ------
%D = 0.04*0.01.*Qmod.*(1-Qmod);
    
iter=1;
%preallocating for speed
AAA = zeros(n*m);
DDD = zeros(n*m,1);
% fprintf('lambda_a(1,1) = %12.10f\n',lambda_a(1,1))
% fprintf('lambda_o(1,1) = %12.10f\n',lambda_o(1,1))
% fprintf('lambda(1,1) = %12.10f\n',lambda(1,1))
% fprintf('f(1,1) = %12.10f\n',f(1,1))
%  pause

while(iter<=(m)*(n-1)+1)
    cnt = (iter -1)/m; %% cnt = 0,1,2,3.... for iter = 1, m+1, 2m+1, 3m+1 etc...
    BB=zeros(n,m);AA=BB;CC=BB;DD=zeros(m,1);
    for i=1:1:m  %% i = x direction discretization formula and columns ie 2nd coordinate in matrix code
        for j=1:1:n %% j = y direction discretization formula and rows ie 1st coordinate in matrix code
            if j==i
                if iter == 1            % 1st or bottom row of grid
                    if i == 1           % 1st or left column
                        DD(i) = (Qmod(cnt+1,i)/dtcal(cnt+1,i)) + G1*(1-f(cnt+1,i))...
                            +((Dg(cnt+1,i)+Dg(cnt+1,i+1))/(dx^2) + (Dg(cnt+2,i)+Dg(cnt+1,i))/(dy^2))*G(cnt+1,i)...
                            -(Dg(cnt+1,i)+Dg(cnt+1,i+1))/(dx^2)*G(cnt+1,i+1)...
                            -(Dg(cnt+1,i)+Dg(cnt+2,i))/(dy^2)*G(cnt+2,i); %src
                        CC(j,i) = (Ds(cnt+1,i)+Ds(cnt+2,i))/(dy^2);
                        BB(j,i) = 1/dtcal(cnt+1,i)-(Ds(cnt+1,i)+Ds(cnt+1,i+1))/(dx^2)...
                            -(Ds(cnt+2,i)+Ds(cnt+1,i))/(dy^2);
                        BB(j,i+1) = (Ds(cnt+1,i)+Ds(cnt+1,i+1))/(dx^2);
                    elseif i == m       % last or rightmost column of grid
                        DD(i) = Qmod(cnt+1,i)/dtcal(cnt+1,i) ...
                            +((Dg(cnt+1,i)+Dg(cnt+1,i-1))/(dx^2) + (Dg(cnt+2,i)+Dg(cnt+1,i))/(dy^2))*G(cnt+1,i)...
                            -(Dg(cnt+1,i)+Dg(cnt+1,i-1))/(dx^2)*G(cnt+1,i-1)...
                            -(Dg(cnt+1,i)+Dg(cnt+2,i))/(dy^2)*G(cnt+2,i);
                        BB(j,i-1)= (Ds(cnt+1,i)+Ds(cnt+1,i-1))/(dx^2);
                        BB(j,i) = 1/dtcal(cnt+1,i)-(Ds(cnt+1,i)+Ds(cnt+1,i-1))/(dx^2)...
                            -(Ds(cnt+2,i)+Ds(cnt+1,i))/(dy^2);
                        CC(j,i) = (Ds(cnt+1,i)+Ds(cnt+2,i))/(dy^2);
                    else
                        DD(i) = Qmod(cnt+1,i)/dtcal(cnt+1,i)...
                            -f_c(cnt+1,i)*(u(cnt+1,i)*(C(cnt+1,i+1)-C(cnt+1,i-1))/(2*dx))...
                            -f_g(cnt+1,i)*(u(cnt+1,i)*(G(cnt+1,i+1)-G(cnt+1,i-1))/(2*dx))...
                            +((Dg(cnt+1,i+1)+Dg(cnt+1,i-1)+2*Dg(cnt+1,i))/(2*dx^2)+(Dg(cnt+1,i+1)+Dg(cnt+1,i))/(dy^2))*G(cnt+1,i)...
                            -(Dg(cnt+1,i+1)+Dg(cnt+1,i))/(2*dx^2)*G(cnt+1,i+1)...
                            -(Dg(cnt+1,i-1)+Dg(cnt+1,i))/(2*dx^2)*G(cnt+1,i-1)...
                            -(Dg(cnt+1,i)+Dg(cnt+2,i))/(dy^2)*G(cnt+2,i);
                        BB(j,i) = 1/dtcal(cnt+1,i)-(Ds(cnt+1,i+1)+Ds(cnt+1,i-1)+2*Ds(cnt+1,i))/(2*dx^2)-(Ds(cnt+2,i)+Ds(cnt+1,i))/(dy^2);
                        BB(j,i-1) = (Ds(cnt+1,i-1)+Ds(cnt+1,i))/(2*dx^2);
                        BB(j,i+1) = (Ds(cnt+1,i+1)+Ds(cnt+1,i))/(2*dx^2);
                        CC(j,i) = (Ds(cnt+1,i)+Ds(cnt+2,i))/(dy^2);
                    end
                elseif iter == (m)*(n-1)+1   %topmost row of grid
                    if i == 1                % leftmost column
                        DD(i) = (Qmod(cnt+1,i)/dtcal(cnt+1,i)) ...
                            +((Dg(cnt+1,i)+Dg(cnt+1,i+1))/(dx^2) + (Dg(cnt,i)+Dg(cnt+1,i))/(dy^2))*G(cnt+1,i)...
                            -(Dg(cnt+1,i)+Dg(cnt+1,i+1))/(dx^2)*G(cnt+1,i+1)...
                            -(Dg(cnt+1,i)+Dg(cnt,i))/(dy^2)*G(cnt,i);
                        AA(j,i) = (Ds(cnt+1,i)+Ds(cnt,i))/(dy^2);
                        BB(j,i) = 1/dtcal(cnt+1,i)-(Ds(cnt+1,i)+Ds(cnt+1,i+1))/(dx^2)...
                            -(Ds(cnt,i)+Ds(cnt+1,i))/(dy^2);
                        BB(j,i+1) = (Ds(cnt+1,i)+Ds(cnt+1,i+1))/(dx^2);
                    elseif i == m            % rightmost column
                        DD(i) = Qmod(cnt+1,i)/dtcal(cnt+1,i) ... %+ G1*(f(cnt+1,i)-lambda_o(cnt+1,i))... % - G1*lambda_a(cnt+1,i)/lambda(cnt+1,i)...
                            +((Dg(cnt+1,i)+Dg(cnt+1,i-1))/(dx^2) + (Dg(cnt,i)+Dg(cnt+1,i))/(dy^2))*G(cnt+1,i)...
                            -(Dg(cnt+1,i)+Dg(cnt+1,i-1))/(dx^2)*G(cnt+1,i-1)...
                            -(Dg(cnt+1,i)+Dg(cnt,i))/(dy^2)*G(cnt,i);
                        BB(j,i-1)= (Ds(cnt+1,i)+Ds(cnt+1,i-1))/(dx^2);
                        BB(j,i) = 1/dtcal(cnt+1,i)-(Ds(cnt+1,i)+Ds(cnt+1,i-1))/(dx^2)...
                            -(Ds(cnt,i)+Ds(cnt+1,i))/(dy^2);
                        AA(j,i) = (Ds(cnt+1,i)+Ds(cnt,i))/(dy^2);
                    else
                        DD(i) = Qmod(cnt+1,i)/dtcal(cnt+1,i)...
                            -f_c(cnt+1,i)*(u(cnt+1,i)*(C(cnt+1,i+1)-C(cnt+1,i-1))/(2*dx))...
                            -f_g(cnt+1,i)*(u(cnt+1,i)*(G(cnt+1,i+1)-G(cnt+1,i-1))/(2*dx))...
                            +((Dg(cnt+1,i+1)+Dg(cnt+1,i-1)+2*Dg(cnt+1,i))/(2*dx^2)+(Dg(cnt+1,i+1)+Dg(cnt+1,i))/(dy^2))*G(cnt+1,i)...
                            -(Dg(cnt+1,i+1)+Dg(cnt+1,i))/(2*dx^2)*G(cnt+1,i+1)...
                            -(Dg(cnt+1,i-1)+Dg(cnt+1,i))/(2*dx^2)*G(cnt+1,i-1)...
                            -(Dg(cnt+1,i)+Dg(cnt,i))/(dy^2)*G(cnt,i);
                        BB(j,i) = 1/dtcal(cnt+1,i)-(Ds(cnt+1,i+1)+Ds(cnt+1,i-1)+2*Ds(cnt+1,i))/(2*dx^2)-(Ds(cnt,i)+Ds(cnt+1,i))/(dy^2);
                        BB(j,i-1) = (Ds(cnt+1,i-1)+Ds(cnt+1,i))/(2*dx^2);
                        BB(j,i+1) = (Ds(cnt+1,i+1)+Ds(cnt+1,i))/(2*dx^2);
                        AA(j,i) = (Ds(cnt+1,i)+Ds(cnt,i))/(dy^2);
                    end
                else                         % interior rows of grid
                    if i == 1                % left most column of grid
                        DD(i) = (Qmod(cnt+1,i)/dtcal(cnt+1,i)) ...
                            -f_c(cnt+1,i)*(v(cnt+1,i)*(C(cnt+2,i)-C(cnt,i))/(2*dy))...
                            -f_g(cnt+1,i)*(v(cnt+1,i)*(G(cnt+2,i)-G(cnt,i))/(2*dy))...
                            +((Dg(cnt+1,i)+Dg(cnt+1,i+1))/(dx^2) + (Dg(cnt,i)+2*Dg(cnt+1,i)+Dg(cnt+2,i))/(2*dy^2))*G(cnt+1,i)...
                            -(Dg(cnt+1,i)+Dg(cnt+1,i+1))/(dx^2)*G(cnt+1,i+1)...
                            -(Dg(cnt+1,i)+Dg(cnt+2,i))/(2*dy^2)*G(cnt+2,i)...
                            -(Dg(cnt+1,i)+Dg(cnt,i))/(2*dy^2)*G(cnt,i);
                        AA(j,i) = (Ds(cnt+1,i)+Ds(cnt,i))/(2*dy^2);
                        BB(j,i) = 1/dtcal(cnt+1,i)-(Ds(cnt+1,i)+Ds(cnt+1,i+1))/(dx^2)...
                            -(Ds(cnt,i)+2*Ds(cnt+1,i)+Ds(cnt+2,i))/(2*dy^2);
                        BB(j,i+1) = (Ds(cnt+1,i)+Ds(cnt+1,i+1))/(dx^2);
                        CC(j,i) = (Ds(cnt+1,i)+Ds(cnt+2,i))/(2*dy^2);
                    elseif i == m            % right most column of grid
                        DD(i) = (1/dtcal(cnt+1,i))*Qmod(cnt+1,i) ...
                            -f_c(cnt+1,i)*(v(cnt+1,i)*(C(cnt+2,i)-C(cnt,i))/(2*dy))...
                            -f_g(cnt+1,i)*(v(cnt+1,i)*(G(cnt+2,i)-G(cnt,i))/(2*dy))...
                            +((Dg(cnt+1,i)+Dg(cnt+1,i-1))/(dx^2) + (Dg(cnt,i)+2*Dg(cnt+1,i)+Dg(cnt+2,i))/(2*dy^2))*G(cnt+1,i)...
                            -(Dg(cnt+1,i)+Dg(cnt+1,i-1))/(dx^2)*G(cnt+1,i-1)...
                            -(Dg(cnt+1,i)+Dg(cnt+2,i))/(2*dy^2)*G(cnt+2,i)...
                            -(Dg(cnt+1,i)+Dg(cnt,i))/(2*dy^2)*G(cnt,i);
                        AA(j,i) = (Ds(cnt+1,i)+Ds(cnt,i))/(2*dy^2);
                        BB(j,i) = 1/dtcal(cnt+1,i)-(Ds(cnt+1,i)+Ds(cnt+1,i-1))/(dx^2)...
                            -(Ds(cnt,i)+2*Ds(cnt+1,i)+Ds(cnt+2,i))/(2*dy^2);
                        BB(j,i-1) = (Ds(cnt+1,i)+Ds(cnt+1,i-1))/(dx^2);
                        CC(j,i) = (Ds(cnt+1,i)+Ds(cnt+2,i))/(2*dy^2);
                    else
                        DD(i) =(Qmod(cnt+1,i)/dtcal(cnt+1,i)) ...
                            - f_c(cnt+1,i)*(u(cnt+1,i)*(C(cnt+1,i+1)-C(cnt+1,i-1))/(2*dx)+v(cnt+1,i)*(C(cnt+2,i)-C(cnt,i))/(2*dy))...
                            - f_g(cnt+1,i)*(u(cnt+1,i)*(G(cnt+1,i+1)-G(cnt+1,i-1))/(2*dx)+v(cnt+1,i)*(G(cnt+2,i)-G(cnt,i))/(2*dy))...
                            - (Dg(cnt+1,i+1)/(2*dx*dx)*(G(cnt+1,i+1)-G(cnt+1,i))... 
                            + Dg(cnt+1,i-1)/(2*dx*dx)*(G(cnt+1,i-1)-G(cnt+1,i))... 
                            + Dg(cnt+1,i)/(2*dx*dx)*(G(cnt+1,i-1)+G(cnt+1,i+1)-2*G(cnt+1,i))...
                            + Dg(cnt+2,i)/(2*dy*dy)*(G(cnt+2,i)-G(cnt+1,i))...
                            + Dg(cnt,i)/(2*dy*dy)*(G(cnt,i)-G(cnt+1,i))...
                            + Dg(cnt+1,i)/(2*dy*dy)*(G(cnt+2,i)+G(cnt,i)-2*G(cnt+1,i)));
                        AA(j,i)   =  (Ds(cnt+1-1,i)+Ds(cnt+1,i))/(2*dy^2);
                        CC(j,i)   =  (Ds(cnt+1,i)+Ds(cnt+1+1,i))/(2*dy^2);
                        BB(j,i)   = 1/dtcal(cnt+1,i)-((1/(2*dx^2))*(Ds(cnt+1,i-1)+2*Ds(cnt+1,i)+Ds(cnt+1,i+1))+...
                            (1/(2*dy^2))*(Ds(cnt+1-1,i)+2*Ds(cnt+1,i)+Ds(cnt+1+1,i)));
                        BB(j,i+1) =  (Ds(cnt+1,i)+Ds(cnt+1,i+1))/(2*dx^2);
                        BB(j,i-1) =  (Ds(cnt+1,i-1)+Ds(cnt+1,i))/(2*dx^2);
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
% disp(reshape(DDD,16,16));
% disp(isnan(DDD)');
% QQQ = AAA\DDD;
% Qnew = reshape(QQQ,m,n)';
disp('Entering bicgstab for saturation calculation');
Qnew=(reshape(bicgstab(AAA,DDD,1*10^(-10),400,[],[]),m,n))';


Qnew(Qnew>1) = 1;


% fprintf('Qnew(1,1) = %12.10f\n',Qnew(1,1))
% fprintf('Qnew(n,m) = %12.10f\n',Qnew(n,m))

%Smax=max(max(Qnew)); %disp(Smax);
%Smin=min(min(Qnew)); %disp(Smin);
   
%pause(3)

%% Solving for concentration of polymer
% recompute characteristics for concentration equation and discretize using
% finite difference on recomputed coordinates

%new redefined coordinates for concentration equation
[xmod2,ymod2] = eval_Xsurf_neumann(x,y,Q,Qnew,G,f,f_s,D,pc_s,pc_g,u,v,dtcal,para,2);

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
                if iter == 1    %lowermost row of grid
                    if i == 1   % leftmost point (source) 
                        DD(i)   = G2/Qnew(cnt+1,i)+ Cmod(cnt+1,i)/dtcal(cnt+1,i); %*src/Qnew(cnt+1,i);
                        BB(j,i) = 1/dtcal(cnt+1,i)+ G1/Qnew(cnt+1,i);%1/dtcal(cnt+1,i)+1/Qnew(cnt+1,i);
                    else 
                        DD(i)   = Cmod(cnt+1,i)/dtcal(cnt+1,i);
                        BB(j,i) =  1/dtcal(cnt+1,i);
                    end             
                elseif iter == (m)*(n-1)+1       %topmost row of grid         
                    if i == m     % rightmost point (sink) 
                        DD(i)  = Cmod(cnt+1,i)/dtcal(cnt+1,i); %- G2*f(cnt+1,i)/Qnew(cnt+1,i)
                        BB(j,i)= 1/dtcal(cnt+1,i); %- G1*f(cnt+1,i)/Qnew(cnt+1,i)
                    else
                        DD(i)   = Cmod(cnt+1,i)/dtcal(cnt+1,i);
                        BB(j,i) =  1/dtcal(cnt+1,i);
                    end                                           
                else %interior rows of grid
                    DD(i) = Cmod(cnt+1,i)/dtcal(cnt+1,i);
                    BB(j,i) =  1/dtcal(cnt+1,i);
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

% 
%fprintf('Cnew(1,1) = %12.10f\n',Cnew(1,1))
%fprintf('Cnew(n,m) = %12.10f\n',Cnew(n,m))
%pause(2)
%Cnew(Cnew<0) = 0;

%Cmax=max(max(Cnew)); disp(Cmax);

%% Solving for Surfactant concentration
% recompute characteristics for surfactant equation and discretize using
% finite difference on recomputed coordinates
     
%new redefined coordinates for concentration equation
[xmod3,ymod3] =  eval_Xsurf_neumann(x,y,Q,Qnew,G,f,f_s,D,pc_s,pc_g,u,v,dtcal,para,3);

%bilinear interpolant for sur conc on redefined coordinates
Gmod = interp2(x,y,G,xmod3,ymod3);

%updating coefficients with interpolated Surf conc 
sigmamod = 10.001./(Gmod+1)-0.001;
sigma_g_mod = -10.001./(Gmod+1).^2;
[swr,sor] = compres(sigmamod,u,v,miua);
lambda_a = compmob(Qmod,C,sor,swr,1);
lambda_o = compmob(Qmod,C,sor,swr,0);
lambda = lambda_a + lambda_o;
f = lambda_a./lambda;
D = KK.*lambda_o.*f;

%pause
if g0~=0
    nso = (Qmod-swr)/(1-swr-sor);
    nso_g = (swr_g.*(Qmod+sor-1)+sor_g.*(Qmod-swr))/(1-swr-sor)^2;
    pc = (omega2.*sigma.*sqrt(phi./KK))./((1-nso)^(1/omega1));
    pc_s = pc./(omega1*(1-nso));
    pc_g = (pc./sigmamod).*sigma_g_mod + pc_s.*nso_g; % == 0 when no surfactant
else
    nsw = (Qmod-swr)/(1-swr-sor);
    pc = ((sigma*omega2*sqrt(phi))./(sqrt(KK))./(nsw)^(1/omega1));
    pc_s = - pc./(omega1*nsw); % $$computing \frac{\partial p_c}{\partial s}$$
    pc_g = (pc./sigmamod).*sigma_g_mod + pc_s.*nsw_g; % $$computing \frac{\partial p_c}{\partial \Gamma}$$
end


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
                if iter == 1            % Lower boundary of grid
                    if i == 1           % Bottom left point
                        DD(i)     =  G3/Qnew(cnt+1,i) + Gmod(cnt+1,i)/dtcal(cnt+1,i);%g0
                        CC(j,i)   =  2*F(cnt+1,i)/(dy*dy);
                        BB(j,i)   =  1/dtcal(cnt+1,i) - (2/dx^2+2/dy^2)*F(cnt+1,i) + G1/Qnew(cnt+1,i);%1/dtcal(cnt+1,i) + 1/Qnew(cnt+1,i);
                        BB(j,i+1) =  2*F(cnt+1,i)/(dx*dx);
                    elseif i == m       % Bottom right point
                        DD(i) = Gmod(cnt+1,i)/dtcal(cnt+1,i);
                        CC(j,i)   =  2*F(cnt+1,i)/(dy*dy);
                        BB(j,i)   =  1/dtcal(cnt+1,i) - (2/dx^2+2/dy^2)*F(cnt+1,i);
                        BB(j,i-1) =  2*F(cnt+1,i)/(dx*dx);
                    else                % Lowermost row of grid
                        DD(i) = Gmod(cnt+1,i)/dtcal(cnt+1,i);
                        CC(j,i)   =  2*F(cnt+1,i)/(dy*dy);
                        BB(j,i-1) =  F(cnt+1,i)/(dx*dx);
                        BB(j,i)   =  1/dtcal(cnt+1,i) - (2/dx^2+2/dy^2)*F(cnt+1,i);
                        BB(j,i+1) =  F(cnt+1,i)/(dx*dx);
                    end
                elseif iter == (m)*(n-1)+1  % Upper boundary of grid
                    if i == 1               % Top left point
                        DD(i) = Gmod(cnt+1,i)/dtcal(cnt+1,i);
                        AA(j,i)   =  2*F(cnt+1,i)/(dy*dy);
                        BB(j,i)   =  1/dtcal(cnt+1,i) - (2/dx^2+2/dy^2)*F(cnt+1,i);
                        BB(j,i+1) =  2*F(cnt+1,i)/(dx*dx);
                    elseif i == m           % Top right point
                        DD(i) = Gmod(cnt+1,i)/dtcal(cnt+1,i);
                        AA(j,i)   =  2*F(cnt+1,i)/(dy*dy);
                        BB(j,i)   =  1/dtcal(cnt+1,i) - (2/dx^2+2/dy^2)*F(cnt+1,i);% ...
                                    %-G1/Qnew(cnt+1,i); % ...
                                    %+(G3*lambda_a(cnt+1,i)/lambda(cnt+1,i))/(Qnew(cnt+1,i)*g0);
                        BB(j,i-1) =  2*F(cnt+1,i)/(dx*dx);
                    else                    % Upper most row of grid
                        DD(i) = Gmod(cnt+1,i)/dtcal(cnt+1,i);
                        AA(j,i)   =  2*F(cnt+1,i)/(dy*dy);
                        BB(j,i+1) =  F(cnt+1,i)/(dx*dx);
                        BB(j,i)   =  1/dtcal(cnt+1,i) - (2/dx^2+2/dy^2)*F(cnt+1,i);
                        BB(j,i-1) =  F(cnt+1,i)/(dx*dx);
                    end
                else                        % Interior rows of grid
                    if i == 1               %
                        DD(i) = Gmod(cnt+1,i)/dtcal(cnt+1,i);
                        AA(j,i) = F(cnt+1,i)/(dy*dy);
                        BB(j,i)   =  1/dtcal(cnt+1,i) - (2/dx^2+2/dy^2)*F(cnt+1,i);
                        BB(j,i+1) =  2*F(cnt+1,i)/(dx*dx);
                        CC(j,i) = F(cnt+1,i)/(dy*dy);
                    elseif i == m
                        DD(i) = Gmod(cnt+1,i)/dtcal(cnt+1,i);
                        AA(j,i) = F(cnt+1,i)/(dy*dy);
                        BB(j,i)   =  1/dtcal(cnt+1,i) - (2/dx^2+2/dy^2)*F(cnt+1,i);
                        BB(j,i-1) =  2*F(cnt+1,i)/(dx*dx);
                        CC(j,i) = F(cnt+1,i)/(dy*dy);
                    else
                        DD(i) = Gmod(cnt+1,i)/dtcal(cnt+1,i);
                        AA(j,i)   =  F(cnt+1,i)/(dy*dy);
                        CC(j,i)   =  F(cnt+1,i)/(dy*dy);
                        BB(j,i)   =  1/dtcal(cnt+1,i) - (2/dx^2+2/dy^2)*F(cnt+1,i);
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


Gmin =min(min(Gnew)); disp(Gmin);
Gmax =max((max(Gnew))); disp(Gmax);
Gnew(Gnew>g0) = g0;
Gnew(Gnew<1e-20) = 0;
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
% Volume of oil recovered at production well
ocut = lambda_o(n,m)*src./lambda(n,m); 
% Volume of water produced at production well
wcut = lambda_a(n,m)*src./lambda(n,m); %lambda_a(n,m)*src*Qnew(n,m)./lambda(n,m);
% Oil still in place as a percentage of volume fraction
ROIP = 100*sum(sum(1 - Qnew))/sum(ones(n*m,1));

