function [UU,CC] = nmmoc3(u,v,~,~,S,C,miuo,miuw,para,dt)
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
    
    % debugging
    % xmod(:,1:5);
    % ymod(:,1:5);
    % pause
    
    % bilinear interpolant for saturation on redefined coordinates
    Qmod = interp2(x,y,Q,xmod,ymod);
    % Qmod(:,1:5);
    % b_s = eval_bs(Qmod,miuw,miuo);
    % % optimally redefined coordinates
    % [xmod,ymod] = eval_X(x,y,b_s,umod,vmod,dt);
    % % saturation values on redefined coordinates interpolated from knot values
    % Qmod = interp2(x,y,Q,xmod,ymod,'linear',0);
    
    %%%% Partial derivative of flux function f(s,c) wrt c
    f_c = - (Q.^3.*(1-Q).^3)./((Q.^3+(0.5+C).*(1-Q).^3).^2);
    % [cx,cy] = gradient(C,para.box.dx);
    D = 0.04*0.01.*Qmod.*(1-Qmod);
    Qnew=zeros(n-2,m-2);  % initialize Qnew for speed of computation
    
    %%%% Solve saturation equation -----------
    for i=2:1:m-1  %% i = x direction discretization formula and columns ie 2nd coordinate in matrix code
        for j=2:1:n-1 %% j = y direction discretization formula and rows ie 1st coordinate in matrix code
            
%             Qnew(j-1,i-1) = (1/(1/dtmod(j,i)+(1/(2*dx*dx))*(D(j,i-1)+2*D(j,i)+D(j,i+1))...
%                 +(1/(2*dy*dy))*(D(j-1,i)+2*D(j,i)+D(j+1,i))))*((1/dtmod(j,i))*Qmod(j,i)...
%                 +((D(j,i)+D(j,i-1))/(2*dx*dx))*Q(j,i-1)+((D(j,i)+D(j-1,i))/(2*dy*dy))...
%                 *Q(j-1,i)+((D(j,i)+D(j,i+1))/(2*dx*dx))*Q(j,i+1)+((D(j,i)+D(j+1,i))...
%                 /(2*dy*dy))*Q(j+1,i)-f_c(j,i)*(u(j,i)*((C(j,i+1)-C(j,i-1))/(2*dx))...
%                 +v(j,i)*((C(j+1,i)-C(j-1,i))/(2*dy))));

            Qnew(j-1,i-1) = (dtmod(j,i))*(1/dtmod(j,i)*Qmod(j,i)+((1/(2*dx*dx))*(D(j,i-1)+2*D(j,i)+D(j,i+1))...
                +(1/(2*dy*dy))*(D(j-1,i)+2*D(j,i)+D(j+1,i)))*Q(j,i)...
                +((D(j,i)+D(j,i-1))/(2*dx*dx))*Q(j,i-1)+((D(j,i)+D(j-1,i))/(2*dy*dy))...
                *Q(j-1,i)+((D(j,i)+D(j,i+1))/(2*dx*dx))*Q(j,i+1)+((D(j,i)+D(j+1,i))...
                /(2*dy*dy))*Q(j+1,i)-f_c(j,i)*(u(j,i)*((C(j,i+1)-C(j,i-1))/(2*dx))...
                +v(j,i)*((C(j+1,i)-C(j-1,i))/(2*dy))));
        end
    end
    a1 = Qnew(1,m-2); %right bottom vertex of Q
    b1 = Qnew(n-2,1); %left top vertex of Q
    c1 = Qnew(n-2,m-2); %right top vertex of Q
    Qcal = [0.79 Qnew(1,:) a1;
        Qnew(:,1) Qnew Qnew(:,m-2);
        b1 Qnew(n-2,:) c1];
    Qnew=Qcal;
    
    
    %%%%%%%%%%debugging----------------------------------
%     [j1,i1,v1] = find(Qnew<0);
%     v1'
%     Qnew(j1,i1) = 1/4 * (Qnew(j1-1,i1)+Qnew(j1+1,i1)+Qnew(j1,i1-1)+Qnew(j1,i1+1));
%     Qnew
%     pause
    %%%%%%%%%%%%%----------------------------------------

    %%%redefining Diffusion coefficient with new saturation values
    D = 0.04*0.01.*Qnew.*(1-Qnew);
    %%% Defining g(s,c) = f(s,c)/s
    g = Q.^2./(Qnew.^3+(0.5+C).*(1-Qnew).^3);
    
    %new redefined coordinates for concentration equation
    [xmod2,ymod2,dtmod2] = eval_X(x,y,g,u,v,dtmod);
    
    %%%%%%%%%%%debugging --------------------
    % i=1:1:m; j=1:1:n; flag = zeros(n,m);
    % if dtmod2(j,i) == dtmod(j,i)
    %    flag(j,i) = 1;
    % end
    % flag
    %%%%%%%%%%%%-----------------------------
    
    % bilinear interpolant for saturation on redefined coordinates
    Cmod = interp2(x,y,C,xmod2,ymod2);
    Cnew=zeros(n-2,m-2);  % initialize Cnew for speed of computation
    %%%%% Solving concentration equation -------
    for i=2:1:m-1  %% i = x direction discretization formula and columns ie 2nd coordinate in matrix code
        for j=2:1:n-1 %% j = y direction discretization formula and rows ie 1st coordinate in matrix code
            
            Cnew(j-1,i-1) = Cmod(j,i)+ dtmod2(j,i)*(D(j,i)/Q(j,i))...
                *((1/(4*dx*dx))*(Qnew(j,i+1)-Qnew(j,i-1))*(C(j,i+1)-C(j,i-1))...
                + (1/(4*dy*dy))*(Qnew(j+1,i)-Qnew(j-1,i))*(C(j+1,i)-C(j-1,i)));
            
        end
    end
    
    disp('mmoc calculation done');
    
    a2 = Cnew(1,m-2); %right bottom vertex of Q
    b2 = Cnew(n-2,1); %left top vertex of Q
    c2 = Cnew(n-2,m-2); %right top vertex of Q
    Ccal = [0.05 Cnew(1,:) a2;
        Cnew(:,1) Cnew Cnew(:,m-2);
        b2 Cnew(n-2,:) c2];
    Cnew=Ccal;
%     [j2,i2,v2] = find(Cnew<0);
%     Cnew(j1,i1) = 1/4 * (Cnew(j1-1,i1)+Cnew(j1+1,i1)+Cnew(j1,i1-1)+Cnew(j1,i1+1));
%     Cnew
%     pause
%     not_done = norm(norm(abs(Qnew-Q),inf),inf) > 10^(-6)
%     Q=Qnew
%     C=Cnew
%     pause
% end


% pause
% Qnew=S;
% Qnew(2:n-1,2:m-1) = Q;

% figure(100)
% plot(x(1,:),S(2,:),'Color','blue');
% hold on
% axis([0 1 -0.5 0.9])
% plot(x(1,:),Qnew(2,:)-0.4,'Color','red');


% figure(101)
% plot(x(1,:),S(floor(n-2),:),'Color','blue');
% hold on
% axis([0 1 -0.5 0.9])
% plot(x(1,:),Qnew(floor(n-2),:)-0.4,'Color','red');
UU=Qnew;
CC=Cnew;

