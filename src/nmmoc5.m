function [UU,CC] = nmmoc5(u,v,S,C,miuo,miuw,para,dt,src,s0,c0)
% code to compute solution of saturation equations by Modified Method of
% Characteristics using explicit formulation (Yuan Yi-Rang 1993)
format long
m=size(S,2); n=size(S,1); % Notice m = para.box.m+1 and n = para.box.n+1
dx = para.box.dx; dy = para.box.dy;
Q = S;
%make an array of dt
dtcal = dt*ones(n,m);

% fixed coordinates
[x,y] = meshgrid(para.box.left:para.box.dx:para.box.right,...
    para.box.bottom:para.box.dy:para.box.top);

% defining coefficient of hyperbolic part, b(s)
%%%% b(s) = df/ds(s) for flow without polymer
% b_s = eval_bs(Q,C,miuw,miuo,0);

% not_done =true;
% while not_done
    
%%%% b(s) = (del/del s)f(s,c) = (s^3)/(s^3+(0.5+c)*(1-s)^3) for polymer flow
b_s = eval_bs(Q,C,miuw,miuo,1);

%determine the correct dt so that backward differencing along
%characteristic tangent does not fall outside domain

% redefined coordinates
[xmod,ymod,dtmod] = eval_X(x,y,b_s,u,v,dtcal);


% bilinear interpolant for saturation on redefined coordinates
Qmod = interp2(x,y,Q,xmod,ymod);

% b_s = eval_bs(Qmod,miuw,miuo);
% % optimally redefined coordinates
% [xmod,ymod] = eval_X(x,y,b_s,umod,vmod,dt);
% % saturation values on redefined coordinates interpolated from knot values
% Qmod = interp2(x,y,Q,xmod,ymod,'linear',0);

%%%% Partial derivative of flux function f(s,c) wrt c
f_c = - (Q.^3.*(1-Q).^3)./((Q.^3+(0.5+C).*(1-Q).^3).^2);
% [cx,cy] = gradient(C,para.box.dx);

%%%%%%%% Mobilities Lambda = (k_relative/viscosity)
%%%%Trial-----
lambda_a = Qmod.^3./(miuo*(0.5+C));
%     lambda_a = Q(:,:,1).^3./(miuw*(1+alpha*(Q(:,:,2)./Q(:,:,1))));
lambda_o = (1-Qmod).^3./miuo;
%%%%%%% -------

lambda = lambda_a +lambda_o;


D = 0.04*0.01.*Qmod.*(1-Qmod);
    
iter=1;
DDD=[];

while(iter<=(m)*(n-1)+1)
    cnt = (iter -1)/m; %% cnt = 0,1,2,3.... for iter = 1, m+1, 2m+1, 3m+1 etc...
    BB=zeros(n,m);AA=BB;CC=BB;DD=zeros(m,1);
    for i=1:1:m  %% i = x direction discretization formula and columns ie 2nd coordinate in matrix code
        for j=1:1:n %% j = y direction discretization formula and rows ie 1st coordinate in matrix code
            if j==i
                if iter == 1
                    if i == 1
                        DD(i) = 1;
                        BB(j,i) =1;
%                         DD(i) = (1/dt)*Qmod(cnt+1,i)-f_c(cnt+1,i)*...
%                             (u(cnt+1,i)*(C(cnt+1,i+1)-C(cnt+1,i))/(dx)+v(cnt+1,i)*(C(cnt+2,i)-C(cnt+1,i))/(dy));%+(1-lambda_a(cnt+1,i)/lambda(cnt+1,i))*src
%                         BB(j,i)   =  1/dt+((1/(dx*dx))*(D(cnt+1,i)-D(cnt+1,i+1))+(1/(dy*dy))*(D(cnt+1,i)-D(cnt+2,i)));
%                         BB(j,i+1) =  (D(cnt+1,i+1)-D(cnt+1,i))/(dx*dx);
%                         CC(j,i)   =  (D(cnt+2,i)-D(cnt+1,i))/(dy*dy);
                    elseif i == m
                        DD(i) = Qmod(cnt+1,i-1)/2+Qmod(cnt+2,i)/2;
                        BB(j,i) = 1;
%                         DD(i) = (1/dt)*Qmod(cnt+1,i)-f_c(cnt+1,i)*(u(cnt+1,i)*(C(cnt+1,i)-C(cnt+1,i-1))/(dx)+...
%                             v(cnt+1,i)*(C(cnt+2,i)-C(cnt+1,i))/(dy));
%                         BB(j,i-1) =  (D(cnt+1,i-1)-D(cnt+1,i))/(dx*dx);
%                         BB(j,i)   =  1/dt+((1/(dx*dx))*(D(cnt+1,i)-D(cnt+1,i-1))+(1/(dy*dy))*(D(cnt+1,i)-D(cnt+2,i)));
%                         CC(j,i)   =  (D(cnt+2,i)-D(cnt+1,i))/(dy*dy);
                    else
                        DD(i) = Qmod(cnt+2,i);
                        BB(j,i) = 1;
%                         DD(i) = (1/dt)*Qmod(cnt+1,i)-f_c(cnt+1,i)*(u(cnt+1,i)*(C(cnt+1,i+1)-C(cnt+1,i-1))/(2*dx)+...
%                             v(cnt+1,i)*(C(cnt+2,i)-C(cnt+1,i))/(dy));
%                         BB(j,i)   = 1/dt-((1/(2*dx*dx))*(D(cnt+1,i-1)+D(cnt+1,i+1))+(1/(dy*dy))*(D(cnt+2,i)));
%                         BB(j,i-1) =  (D(cnt+1,i-1)+D(cnt+1,i))/(2*dx*dx);
%                         BB(j,i+1) =  (D(cnt+1,i)+D(cnt+1,i+1))/(2*dx*dx);
%                         CC(j,i)   =  (D(cnt+2,i)-D(cnt+1,i))/(dy*dy);
                    end
                elseif iter == (m)*(n-1)+1
                    if i == 1
                        DD(i) = Qmod(cnt+1,i+1)/2 + Qmod(cnt,i)/2;
                        BB(j,i)= 1;
%                         DD(i) = (1/dt)*Qmod(cnt+1,i)-f_c(cnt+1,i)*(u(cnt+1,i)*(C(cnt+1,i+1)-C(cnt+1,i))/(dx)+...
%                             v(cnt+1,i)*(C(cnt+1,i)-C(cnt,i))/(dy));
%                         AA(j,i)   =  (D(cnt,i)-D(cnt+1,i))/(dy*dy);
%                         BB(j,i)   =  1/dt+((1/(dx*dx))*(D(cnt+1,i)-D(cnt+1,i+1))+(1/(dy*dy))*(D(cnt+1,i)-D(cnt,i)));
%                         BB(j,i+1) =  (D(cnt+1,i+1)-D(cnt+1,i))/(dx*dx);
                    elseif i == m
                        DD(i) = Qmod(cnt+1,i-1)/2 + Qmod(cnt,i)/2;
                        BB(j,i) = 1;
%                         DD(i) = (1/dt)*Qmod(cnt+1,i)-f_c(cnt+1,i)*(u(cnt+1,i)*(C(cnt+1,i)-C(cnt+1,i-1))/(dx)+...
%                             v(cnt+1,i)*(C(cnt+1,i)-C(cnt,i))/(dy));
%                         AA(j,i)   =  (D(cnt,i)-D(cnt+1,i))/(dy*dy);
%                         BB(j,i-1) =  (D(cnt+1,i-1)-D(cnt+1,i))/(dx*dx);
%                         BB(j,i)   = 1/dt-((1/(dx*dx))*(D(cnt+1,i-1)-D(cnt+1,i))+(1/(dy*dy))*(D(cnt,i)-D(cnt+1,i)));
                    else
                        DD(i) = Qmod(cnt,i);
                        BB(j,i) = 1;
%                         DD(i) = (1/dt)*Qmod(cnt+1,i)-f_c(cnt+1,i)*(u(cnt+1,i)*(C(cnt+1,i+1)-C(cnt+1,i-1))/(2*dx)+...
%                             v(cnt+1,i)*(C(cnt+1,i)-C(cnt,i))/(dy));
%                         AA(j,i)   =  (D(cnt,i)-D(cnt+1,i))/(dy*dy);
%                         BB(j,i)   = 1/dt-((1/(2*dx*dx))*(D(cnt+1,i-1)+D(cnt+1,i+1))+(1/(dy*dy))*(D(cnt,i)));
%                         BB(j,i-1) =  (D(cnt+1,i-1)+D(cnt+1,i))/(2*dx*dx);
%                         BB(j,i+1) =  (D(cnt+1,i)+D(cnt+1,i+1))/(2*dx*dx);
                    end
                else
                    if i == 1
                        DD(i) = Qmod(cnt+1,i+1);
                        BB(j,i) = 1;
%                         DD(i) = (1/dt)*Qmod(cnt+1,i)-f_c(cnt+1,i)*(u(cnt+1,i)*(C(cnt+1,i+1)-C(cnt+1,i))/(dx)+...
%                             v(cnt+1,i)*(C(cnt+2,i)-C(cnt,i))/(2*dy));
%                         AA(j,i)   =  (D(cnt,i)+D(cnt+1,i))/(2*dy*dy);
%                         BB(j,i)   = 1/dt-((1/(dx*dx))*(D(cnt+1,i+1))+(1/(2*dy*dy))*(D(cnt,i)+D(cnt+2,i)));
%                         BB(j,i+1) =  (D(cnt+1,i+1)-D(cnt+1,i))/(dx*dx);
%                         CC(j,i)   =  (D(cnt+1,i)+D(cnt+2,i))/(2*dy*dy);
                    elseif i == m
                        DD(i) = Qmod(cnt+1,i-1);
                        BB(j,i) = 1;
%                         DD(i) = (1/dt)*Qmod(cnt+1,i)-f_c(cnt+1,i)*(u(cnt+1,i)*(C(cnt+1,i)-C(cnt+1,i-1))/(dx)+...
%                             v(cnt+1,i)*(C(cnt+2,i)-C(cnt,i))/(2*dy));
%                         AA(j,i)   =  (D(cnt,i)+D(cnt+1,i))/(2*dy*dy);
%                         BB(j,i-1) =  (D(cnt+1,i-1)-D(cnt+1,i))/(dx*dx);
%                         BB(j,i)   = 1/dt-((1/(dx*dx))*(D(cnt+1,i-1))+(1/(2*dy*dy))*(D(cnt,i)+D(cnt+2,i)));
%                         CC(j,i)   =  (D(cnt+1,i)+D(cnt+2,i))/(2*dy*dy);
                    else
                        DD(i) = (1/dt)*Qmod(cnt+1,i)-f_c(cnt+1,i)*(u(cnt+1,i)*(C(cnt+1,i+1)-C(cnt+1,i-1))/(2*dx)+...
                            v(cnt+1,i)*(C(cnt+2,i)-C(cnt,i))/(2*dy));
                        AA(j,i)   =  (D(cnt+1-1,i)+D(cnt+1,i))/(2*dy*dy);
                        CC(j,i)   =  (D(cnt+1,i)+D(cnt+1+1,i))/(2*dy*dy);
                        BB(j,i)   = 1/dt-((1/(2*dx*dx))*(D(cnt+1,i-1)+2*D(cnt+1,i)+D(cnt+1,i+1))+...
                            (1/(2*dy*dy))*(D(cnt+1-1,i)+2*D(cnt+1,i)+D(cnt+1+1,i)));
                        BB(j,i+1) =  (D(cnt+1,i)+D(cnt+1,i+1))/(2*dx*dx);
                        BB(j,i-1) =  (D(cnt+1,i-1)+D(cnt+1,i))/(2*dx*dx);
                    end
                end
                                
            end
        end
    end
    if(iter==1)
        AAA = [BB,CC,zeros(n,(n)*(m)-2*(m))];
    elseif(iter==(n-1)*(m)+1)
        AAA = [AAA;zeros(n,(n)*(m)-2*(m)),AA,BB];
    else
        AAA = [AAA;zeros(n,m*(cnt-1)),AA,BB,CC,zeros(n,(n)*(m)-3*(m)-m*(cnt-1))];
    end
    DDD = [DDD;DD];
    iter = iter+m;
end


% QQQ = AAA\DDD;
% Qnew = reshape(QQQ,m,n)';

Qnew=(reshape(bicgstab(AAA,DDD,1*10^(-10),600,[],[]),m,n))';
% Qnew(1,1)=s0;
% disp('mmoc calculation done');


    
    %%%%%%%%%%debugging----------------------------------
%     [j1,i1,v1] = find(Qnew<0);
%     v1'
%     Qnew(j1,i1) = 1/4 * (Qnew(j1-1,i1)+Qnew(j1+1,i1)+Qnew(j1,i1-1)+Qnew(j1,i1+1));
%     Qnew
%     pause
    %%%%%%%%%%%%%----------------------------------------
    
%     lambda_a = Qmod.^3./(miuo*(0.5+C));
%     lambda_o = (1-Qmod).^3./miuo;
%     %%%%%%% -------
%     
%     lambda = lambda_a +lambda_o;



%%%redefining Diffusion coefficient with new saturation values
D = 0.04*0.01.*Qnew.*(1-Qnew);
%%% Defining g(s,c) = f(s,c)/s
g = Qmod.^2./((Qmod.^3+(0.5+C).*(1-Qmod).^3));

%new redefined coordinates for concentration equation
[xmod2,ymod2,dtmod2] = eval_X(x,y,g,u,v,dtmod);

%bilinear interpolant for saturation on redefined coordinates
Cmod = interp2(x,y,C,xmod2,ymod2);
% Cnew=zeros(n,m);  % initialize Cnew for speed of computation

iter =1;
DDD=[];

while(iter<=(m)*(n-1)+1)
    cnt = (iter -1)/m; %% cnt = 0,1,2,3.... for iter = 1, m+1, 2m+1, 3m+1 etc...
    BB=zeros(n,m);AA=BB;CC=BB;DD=zeros(m,1);
    for i=1:1:m  %% i = x direction discretization formula and columns ie 2nd coordinate in matrix code
        for j=1:1:n %% j = y direction discretization formula and rows ie 1st coordinate in matrix code
            if j==i
                if iter == 1
                    if i == 1
                        DD(i) = c0;
                        BB(j,i) =1;
%                         DD(i) = (1/dt)*Qmod(cnt+1,i)-f_c(cnt+1,i)*...
%                             (u(cnt+1,i)*(C(cnt+1,i+1)-C(cnt+1,i))/(dx)+v(cnt+1,i)*(C(cnt+2,i)-C(cnt+1,i))/(dy));%+(1-lambda_a(cnt+1,i)/lambda(cnt+1,i))*src
%                         BB(j,i)   =  1/dt+((1/(dx*dx))*(D(cnt+1,i)-D(cnt+1,i+1))+(1/(dy*dy))*(D(cnt+1,i)-D(cnt+2,i)));
%                         BB(j,i+1) =  (D(cnt+1,i+1)-D(cnt+1,i))/(dx*dx);
%                         CC(j,i)   =  (D(cnt+2,i)-D(cnt+1,i))/(dy*dy);
                    elseif i == m
                        DD(i) = Cmod(cnt+1,i-1)/2+Cmod(cnt+2,i)/2;
                        BB(j,i) = 1;
%                         DD(i) = (1/dt)*Qmod(cnt+1,i)-f_c(cnt+1,i)*(u(cnt+1,i)*(C(cnt+1,i)-C(cnt+1,i-1))/(dx)+...
%                             v(cnt+1,i)*(C(cnt+2,i)-C(cnt+1,i))/(dy));
%                         BB(j,i-1) =  (D(cnt+1,i-1)-D(cnt+1,i))/(dx*dx);
%                         BB(j,i)   =  1/dt+((1/(dx*dx))*(D(cnt+1,i)-D(cnt+1,i-1))+(1/(dy*dy))*(D(cnt+1,i)-D(cnt+2,i)));
%                         CC(j,i)   =  (D(cnt+2,i)-D(cnt+1,i))/(dy*dy);
                    else
                        DD(i) = Cmod(cnt+2,i);
                        BB(j,i) = 1;
%                         DD(i) = (1/dt)*Qmod(cnt+1,i)-f_c(cnt+1,i)*(u(cnt+1,i)*(C(cnt+1,i+1)-C(cnt+1,i-1))/(2*dx)+...
%                             v(cnt+1,i)*(C(cnt+2,i)-C(cnt+1,i))/(dy));
%                         BB(j,i)   = 1/dt-((1/(2*dx*dx))*(D(cnt+1,i-1)+D(cnt+1,i+1))+(1/(dy*dy))*(D(cnt+2,i)));
%                         BB(j,i-1) =  (D(cnt+1,i-1)+D(cnt+1,i))/(2*dx*dx);
%                         BB(j,i+1) =  (D(cnt+1,i)+D(cnt+1,i+1))/(2*dx*dx);
%                         CC(j,i)   =  (D(cnt+2,i)-D(cnt+1,i))/(dy*dy);
                    end
                elseif iter == (m)*(n-1)+1
                    if i == 1
                        DD(i) = Cmod(cnt+1,i+1)/2 + Cmod(cnt,i)/2;
                        BB(j,i)= 1;
%                         DD(i) = (1/dt)*Qmod(cnt+1,i)-f_c(cnt+1,i)*(u(cnt+1,i)*(C(cnt+1,i+1)-C(cnt+1,i))/(dx)+...
%                             v(cnt+1,i)*(C(cnt+1,i)-C(cnt,i))/(dy));
%                         AA(j,i)   =  (D(cnt,i)-D(cnt+1,i))/(dy*dy);
%                         BB(j,i)   =  1/dt+((1/(dx*dx))*(D(cnt+1,i)-D(cnt+1,i+1))+(1/(dy*dy))*(D(cnt+1,i)-D(cnt,i)));
%                         BB(j,i+1) =  (D(cnt+1,i+1)-D(cnt+1,i))/(dx*dx);
                    elseif i == m
                        DD(i) = Cmod(cnt+1,i-1)/2 + Cmod(cnt,i)/2;
                        BB(j,i) = 1;
%                         DD(i) = (1/dt)*Qmod(cnt+1,i)-f_c(cnt+1,i)*(u(cnt+1,i)*(C(cnt+1,i)-C(cnt+1,i-1))/(dx)+...
%                             v(cnt+1,i)*(C(cnt+1,i)-C(cnt,i))/(dy));
%                         AA(j,i)   =  (D(cnt,i)-D(cnt+1,i))/(dy*dy);
%                         BB(j,i-1) =  (D(cnt+1,i-1)-D(cnt+1,i))/(dx*dx);
%                         BB(j,i)   = 1/dt-((1/(dx*dx))*(D(cnt+1,i-1)-D(cnt+1,i))+(1/(dy*dy))*(D(cnt,i)-D(cnt+1,i)));
                    else
                        DD(i) = Cmod(cnt,i);
                        BB(j,i) = 1;
%                         DD(i) = (1/dt)*Qmod(cnt+1,i)-f_c(cnt+1,i)*(u(cnt+1,i)*(C(cnt+1,i+1)-C(cnt+1,i-1))/(2*dx)+...
%                             v(cnt+1,i)*(C(cnt+1,i)-C(cnt,i))/(dy));
%                         AA(j,i)   =  (D(cnt,i)-D(cnt+1,i))/(dy*dy);
%                         BB(j,i)   = 1/dt-((1/(2*dx*dx))*(D(cnt+1,i-1)+D(cnt+1,i+1))+(1/(dy*dy))*(D(cnt,i)));
%                         BB(j,i-1) =  (D(cnt+1,i-1)+D(cnt+1,i))/(2*dx*dx);
%                         BB(j,i+1) =  (D(cnt+1,i)+D(cnt+1,i+1))/(2*dx*dx);
                    end
                else
                    if i == 1
                        DD(i) = Cmod(cnt+1,i+1);
                        BB(j,i) = 1;
%                         DD(i) = (1/dt)*Qmod(cnt+1,i)-f_c(cnt+1,i)*(u(cnt+1,i)*(C(cnt+1,i+1)-C(cnt+1,i))/(dx)+...
%                             v(cnt+1,i)*(C(cnt+2,i)-C(cnt,i))/(2*dy));
%                         AA(j,i)   =  (D(cnt,i)+D(cnt+1,i))/(2*dy*dy);
%                         BB(j,i)   = 1/dt-((1/(dx*dx))*(D(cnt+1,i+1))+(1/(2*dy*dy))*(D(cnt,i)+D(cnt+2,i)));
%                         BB(j,i+1) =  (D(cnt+1,i+1)-D(cnt+1,i))/(dx*dx);
%                         CC(j,i)   =  (D(cnt+1,i)+D(cnt+2,i))/(2*dy*dy);
                    elseif i == m
                        DD(i) = Cmod(cnt+1,i-1);
                        BB(j,i) = 1;
%                         DD(i) = (1/dt)*Qmod(cnt+1,i)-f_c(cnt+1,i)*(u(cnt+1,i)*(C(cnt+1,i)-C(cnt+1,i-1))/(dx)+...
%                             v(cnt+1,i)*(C(cnt+2,i)-C(cnt,i))/(2*dy));
%                         AA(j,i)   =  (D(cnt,i)+D(cnt+1,i))/(2*dy*dy);
%                         BB(j,i-1) =  (D(cnt+1,i-1)-D(cnt+1,i))/(dx*dx);
%                         BB(j,i)   = 1/dt-((1/(dx*dx))*(D(cnt+1,i-1))+(1/(2*dy*dy))*(D(cnt,i)+D(cnt+2,i)));
%                         CC(j,i)   =  (D(cnt+1,i)+D(cnt+2,i))/(2*dy*dy);

%%%%%%%%%%%%%%%%% ACTUAL CODE BELOW %%%%%%%%%%%%%%%%%%%%%%
%                 if iter == 1
%                     if i == 1
%                         DD(i) = (1/dt)*Cmod(cnt+1,i); %-(1-lambda_a(cnt+1,i)/lambda(cnt+1,i))*src/Qnew(cnt+1,i)
%                         BB(j,i)   =  1/dt-(1/(dx*dx))*(D(cnt+1,i)/Qnew(cnt+1,i))*(Qnew(cnt+1,i+1)-Qnew(cnt+1,i))...
%                             -(1/(dy*dy))*(D(cnt+1,i)/Qnew(cnt+1,i))*(Qnew(cnt+2,i)-Qnew(cnt+1,i));
%                         BB(j,i+1) =  (D(cnt+1,i)/Qnew(cnt+1,i))*(Qnew(cnt+1,i+1)-Qnew(cnt+1,i))/(dx*dx);
%                         CC(j,i)   =  (D(cnt+1,i)/Qnew(cnt+1,i))*(Qnew(cnt+2,i)-Qnew(cnt+1,i))/(dy*dy);
%                     elseif i == m
%                         DD(i) = (1/dt)*Cmod(cnt+1,i);
%                         BB(j,i-1) =  (D(cnt+1,i)/Qnew(cnt+1,i))*(Qnew(cnt+1,i-1)-Qnew(cnt+1,i))/(dx*dx);
%                         BB(j,i)   =  1/dt+(D(cnt+1,i)/Qnew(cnt+1,i))*(1/(dx*dx))*(Qnew(cnt+1,i)-Qnew(cnt+1,i-1))...
%                             -(D(cnt+1,i)/Qnew(cnt+1,i))*(1/(dy*dy))*(Qnew(cnt+2,i)-Qnew(cnt+1,i));
%                         CC(j,i)   =  (D(cnt+1,i)/Qnew(cnt+1,i))*(Qnew(cnt+2,i)-Qnew(cnt+1,i))/(dy*dy);
%                     else
%                         DD(i) = (1/dt)*Cmod(cnt+1,i);
%                         BB(j,i)   = 1/dt-(D(cnt+1,i)/Qnew(cnt+1,i))*(1/(dy*dy))*(Qnew(cnt+2,i)-Qnew(cnt+1,i));
%                         BB(j,i-1) =  (D(cnt+1,i)/Qnew(cnt+1,i))*(Qnew(cnt+1,i-1)-Qnew(cnt+1,i+1))/(4*dx*dx);
%                         BB(j,i+1) =  (D(cnt+1,i)/Qnew(cnt+1,i))*(Qnew(cnt+1,i+1)-Qnew(cnt+1,i-1))/(4*dx*dx);
%                         CC(j,i)   =  (D(cnt+1,i)/Qnew(cnt+1,i))*(Qnew(cnt+2,i)-Qnew(cnt+1,i))/(dy*dy);
%                     end
%                 elseif iter == (m)*(n-1)+1
%                     if i == 1
%                         DD(i) = (1/dt)*Cmod(cnt+1,i);
%                         AA(j,i)   =  (D(cnt+1,i)/Qnew(cnt+1,i))*(Qnew(cnt,i)-Qnew(cnt+1,i))/(dy*dy);
%                         BB(j,i)   =  1/dt+(D(cnt+1,i)/Qnew(cnt+1,i))*(1/(dx*dx))*(Qnew(cnt+1,i)-Qnew(cnt+1,i+1))...
%                             +(D(cnt+1,i)/Qnew(cnt+1,i))*(1/(dy*dy))*(Qnew(cnt+1,i)-Qnew(cnt,i));
%                         BB(j,i+1) =  (D(cnt+1,i)/Qnew(cnt+1,i))*(Qnew(cnt+1,i+1)-Qnew(cnt+1,i))/(dx*dx);
%                     elseif i == m
%                         DD(i) = (1/dt)*Cmod(cnt+1,i);
%                         AA(j,i)   =  (D(cnt+1,i)/Qnew(cnt+1,i))*(Qnew(cnt,i)-Qnew(cnt+1,i))/(dy*dy);
%                         BB(j,i-1) =  (D(cnt+1,i)/Qnew(cnt+1,i))*(Qnew(cnt+1,i-1)-Qnew(cnt+1,i))/(dx*dx);
%                         BB(j,i)   = 1/dt+(D(cnt+1,i)/Qnew(cnt+1,i))*(1/(dx*dx))*(Qnew(cnt+1,i)-Qnew(cnt+1,i-1))...
%                             +(D(cnt+1,i)/Qnew(cnt+1,i))*(1/(dy*dy))*(Qnew(cnt+1,i)-Qnew(cnt,i));
%                     else
%                         DD(i) = (1/dt)*Cmod(cnt+1,i);
%                         AA(j,i)   =  (D(cnt+1,i)/Qnew(cnt+1,i))*(Qnew(cnt,i)-Qnew(cnt+1,i))/(dy*dy);
%                         BB(j,i)   = 1/dt+(D(cnt+1,i)/Qnew(cnt+1,i))*(1/(dy*dy))*(Qnew(cnt+1,i)-Qnew(cnt,i));
%                         BB(j,i-1) =  (D(cnt+1,i)/Qnew(cnt+1,i))*(Qnew(cnt+1,i-1)-Qnew(cnt+1,i+1))/(4*dx*dx);
%                         BB(j,i+1) =  (D(cnt+1,i)/Qnew(cnt+1,i))*(Qnew(cnt+1,i+1)-Qnew(cnt+1,i-1))/(4*dx*dx);
%                     end
%                 else
%                     if i == 1
%                         DD(i) = (1/dt)*Cmod(cnt+1,i);
%                         AA(j,i)   =  (D(cnt+1,i)/Qnew(cnt+1,i))*(Qnew(cnt,i)-Qnew(cnt+2,i))/(4*dy*dy);
%                         BB(j,i)   = 1/dt-D(cnt+1,i)/Qnew(cnt+1,i)*(1/(dx*dx))*(Qnew(cnt+1,i+1)-Qnew(cnt+1,i));
%                         BB(j,i+1) =  (D(cnt+1,i)/Qnew(cnt+1,i))*(Qnew(cnt+1,i+1)-Qnew(cnt+1,i))/(dx*dx);
%                         CC(j,i)   =  (D(cnt+1,i)/Qnew(cnt+1,i))*(Qnew(cnt+2,i)-Qnew(cnt,i))/(4*dy*dy);
%                     elseif i == m
%                         DD(i) = (1/dt)*Cmod(cnt+1,i);
%                         AA(j,i)   =  D(cnt+1,i)/Qnew(cnt+1,i)*(Qnew(cnt,i)-Qnew(cnt+2,i))/(4*dy*dy);
%                         BB(j,i-1) =  D(cnt+1,i)/Qnew(cnt+1,i)*(Qnew(cnt+1,i-1)-Qnew(cnt+1,i))/(dx*dx);
%                         BB(j,i)   = 1/dt+D(cnt+1,i)/Qnew(cnt+1,i)*(1/(dx*dx))*(Qnew(cnt+1,i)-Qnew(cnt+1,i-1));
%                         CC(j,i)   =  D(cnt+1,i)/Qnew(cnt+1,i)*(Qnew(cnt+2,i)-Qnew(cnt,i))/(4*dy*dy);
%%%%%%%%%%%%%% ACTUAL CODE ABOVE %%%%%%%%%%%%%%%%%%%%%%%
                    else
                                              
                        DD(i) = (1/dt)*Cmod(cnt+1,i);
                        AA(j,i)   =  -D(cnt+1,i)/Qnew(cnt+1,i)*(Q(cnt+2,i)-Q(cnt,i))/(4*dy*dy);
                        CC(j,i)   =  D(cnt+1,i)/Qnew(cnt+1,i)*(Q(cnt+2,i)-Q(cnt,i))/(4*dy*dy);
                        BB(j,i)   =  1/dt;
                        BB(j,i+1) =  D(cnt+1,i)/Qnew(cnt+1,i)*(Q(cnt+1,i+1)-Q(cnt+1,i-1))/(4*dx*dx);
                        BB(j,i-1) =  -D(cnt+1,i)/Qnew(cnt+1,i)*(Q(cnt+1,i+1)-Q(cnt+1,i-1))/(4*dx*dx);
                    end
                end
                                
            end
        end
    end
    if(iter==1)
        AAA = [BB,CC,zeros(n,(n)*(m)-2*(m))];
    elseif(iter==(n-1)*(m)+1)
        AAA = [AAA;zeros(n,(n)*(m)-2*(m)),AA,BB];
    else
        AAA = [AAA;zeros(n,m*(cnt-1)),AA,BB,CC,zeros(n,(n)*(m)-3*(m)-m*(cnt-1))];
    end
    DDD = [DDD;DD];
    iter = iter+m;
end

Cnew=(reshape(bicgstab(AAA,DDD,1*10^(-10),600,[],[]),m,n))';
% Cnew(1,1)=c0;
disp('mmoc calculation done');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%     %%%redefining Diffusion coefficient with new saturation values
%     D = 0.04*0.01.*Qnew.*(1-Qnew);
%     %%% Defining g(s,c) = f(s,c)/s
%     g = Qnew.^2./(Qnew.^3+(0.5+C).*(1-Qnew).^3);
%     
%     %new redefined coordinates for concentration equation
%     [xmod2,ymod2,dtmod2] = eval_X(x,y,g,u,v,dtmod);
%     

%     
%     % bilinear interpolant for saturation on redefined coordinates
%     Cmod = interp2(x,y,C,xmod2,ymod2);
%     Cnew=zeros(n-2,m-2);  % initialize Cnew for speed of computation
%     %%%%% Solving concentration equation -------
%     for i=2:1:m-1  %% i = x direction discretization formula and columns ie 2nd coordinate in matrix code
%         for j=2:1:n-1 %% j = y direction discretization formula and rows ie 1st coordinate in matrix code
%             
%             Cnew(j-1,i-1) = Cmod(j,i)+ dtmod2(j,i)*(D(j,i)/Qnew(j,i))...
%                 *((1/(4*dx*dx))*(Qnew(j,i+1)-Qnew(j,i-1))*(C(j,i+1)-C(j,i-1))...
%                 + (1/(4*dy*dy))*(Qnew(j+1,i)-Qnew(j-1,i))*(C(j+1,i)-C(j-1,i)));
%             
%         end
%     end
%     
%     disp('mmoc calculation done');
%     
%     a2 = Cnew(1,m-2); %right bottom vertex of Q
%     b2 = Cnew(n-2,1); %left top vertex of Q
%     c2 = Cnew(n-2,m-2); %right top vertex of Q
%     Ccal = [0.05 Cnew(1,:) a2;
%         Cnew(:,1) Cnew Cnew(:,m-2);
%         b2 Cnew(n-2,:) c2];
%     Cnew=Ccal;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %     [j2,i2,v2] = find(Cnew<0);
% %     Cnew(j1,i1) = 1/4 * (Cnew(j1-1,i1)+Cnew(j1+1,i1)+Cnew(j1,i1-1)+Cnew(j1,i1+1));
% %     Cnew
% %     pause
% %     not_done = norm(norm(abs(Qnew-Q),inf),inf) > 10^(-6)
% %     Q=Qnew
% %     C=Cnew
% %     pause
% % end
% 
% 
% % pause
% % Qnew=S;
% % Qnew(2:n-1,2:m-1) = Q;
% 
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
Qnew(1:6,1:5);
Cnew(1:6,1:5);
UU=Qnew;
CC=Cnew;

