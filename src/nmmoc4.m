function [UU,CC] = nmmoc4(u,v,S,C,miuo,miuw,para,dt,src,s0,c0)
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
    % Qmod(:,1:5);
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

    Qcal = [Qmod(1,1) Qmod(1,:) Qmod(1,m);
            Qmod(:,1) Qmod Qmod(:,m);
            Qmod(n,1) Qmod(n,:) Qmod(n,m)];
    Ccal = [C(1,1) C(1,:) C(1,m);
            C(:,1) C C(:,m);
            C(n,1) C(n,:) C(n,m)];

    D = 0.04*0.01.*Qcal.*(1-Qcal);
    Qcal = [Q(1,1) Q(1,:) Q(1,m);
            Q(:,1) Q Q(:,m);
            Q(n,1) Q(n,:) Q(n,m)];
    Qnew=zeros(n,m);  % initialize Qnew for speed of computation
    
    
    %%%% Solve saturation equation -----------
    for i=1:1:m  %% i = x direction discretization formula and columns ie 2nd coordinate in matrix code
        for j=1:1:n %% j = y direction discretization formula and rows ie 1st coordinate in matrix code
            
            if i==1 && j==1
%                 Qnew(j,i) = s0;
                Qnew(j,i) = Qmod(j,i)+(dtmod(j,i))*(-((1/(2*dx*dx))*(D(j+1,i)+2*D(j+1,i+1)+D(j+1,i+2))...
                -(1/(2*dy*dy))*(D(j,i+1)+2*D(j+1,i+1)+D(j+2,i+1)))*Qcal(j+1,i+1)...
                +((D(j+1,i+1)+D(j+1,i))/(2*dx*dx))*Qcal(j+1,i)+((D(j+1,i+1)+D(j,i+1))/(2*dy*dy))...
                *Qcal(j,i+1)+((D(j+1,i+1)+D(j+1,i+2))/(2*dx*dx))*Qcal(j+1,i+2)+((D(j+1,i+1)+D(j+2,i+1))...
                /(2*dy*dy))*Qcal(j+2,i+1)-f_c(j,i)*(u(j,i)*((Ccal(j+1,i+2)-Ccal(j+1,i))/(2*dx))...
                +v(j,i)*((Ccal(j+2,i+1)-Ccal(j,i+1))/(2*dy)))+(1-lambda_a(j,i)/lambda(j,i))*src);
            else

                Qnew(j,i) = Qmod(j,i)+(dtmod(j,i))*(-((1/(2*dx*dx))*(D(j+1,i)+2*D(j+1,i+1)+D(j+1,i+2))...
                -(1/(2*dy*dy))*(D(j,i+1)+2*D(j+1,i+1)+D(j+2,i+1)))*Qcal(j+1,i+1)...
                +((D(j+1,i+1)+D(j+1,i))/(2*dx*dx))*Qcal(j+1,i)+((D(j+1,i+1)+D(j,i+1))/(2*dy*dy))...
                *Qcal(j,i+1)+((D(j+1,i+1)+D(j+1,i+2))/(2*dx*dx))*Qcal(j+1,i+2)+((D(j+1,i+1)+D(j+2,i+1))...
                /(2*dy*dy))*Qcal(j+2,i+1)-f_c(j,i)*(u(j,i)*((Ccal(j+1,i+2)-Ccal(j+1,i))/(2*dx))...
                +v(j,i)*((Ccal(j+2,i+1)-Ccal(j,i+1))/(2*dy))));
            end
        end
    end
%     Qnew(1,1)=s0;
    
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


    Ccal = [C(1,1) C(1,:) C(1,m);
            C(:,1) C C(:,m);
            C(n,1) C(n,:) C(n,m)];
        
    %%%redefining Diffusion coefficient with new saturation values
    D = 0.04*0.01.*Qnew.*(1-Qnew);
    %%% Defining g(s,c) = f(s,c)/s
    g = Qnew.^2./(Qnew.^3+(0.5+C).*(1-Qnew).^3);
    
    %new redefined coordinates for concentration equation
    [xmod2,ymod2,dtmod2] = eval_X(x,y,g,u,v,dtmod);
    
    %bilinear interpolant for saturation on redefined coordinates
    Cmod = interp2(x,y,C,xmod2,ymod2);
    Cnew=zeros(n,m);  % initialize Cnew for speed of computation
    %%%%% Solving concentration equation -------
    for i=1:1:m  %% i = x direction discretization formula and columns ie 2nd coordinate in matrix code
        for j=1:1:n %% j = y direction discretization formula and rows ie 1st coordinate in matrix code
            if i==1 && j==1
                Cnew(j,i) = c0;
%                 Cnew(j,i) = Cmod(j,i)+ dtmod2(j,i)*(D(j,i)/Qnew(j,i))...
%                 *((1/(4*dx*dx))*(Qcal(j+1,i+2)-Qcal(j+1,i))*(Ccal(j+1,i+2)-Ccal(j+1,i))...
%                 + (1/(4*dy*dy))*(Qcal(j+2,i+1)-Qcal(j,i+1))*(Ccal(j+2,i+1)-Ccal(j,i+1))+(1-lambda_a(j,i)/lambda(j,i))*src);
            else    
                Cnew(j,i) = Cmod(j,i)+ dtmod2(j,i)*(D(j,i)/Qnew(j,i))...
                    *((1/(4*dx*dx))*(Qcal(j+1,i+2)-Qcal(j+1,i))*(Ccal(j+1,i+2)-Ccal(j+1,i))...
                    + (1/(4*dy*dy))*(Qcal(j+2,i+1)-Qcal(j,i+1))*(Ccal(j+2,i+1)-Ccal(j,i+1)));
            end
        end
    end
    Cnew(1,1)=c0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%     %%%redefining Diffusion coefficient with new saturation values
%     D = 0.04*0.01.*Qnew.*(1-Qnew);
%     %%% Defining g(s,c) = f(s,c)/s
%     g = Qnew.^2./(Qnew.^3+(0.5+C).*(1-Qnew).^3);
%     
%     %new redefined coordinates for concentration equation
%     [xmod2,ymod2,dtmod2] = eval_X(x,y,g,u,v,dtmod);
%     
%     %%%%%%%%%%%debugging --------------------
%     % i=1:1:m; j=1:1:n; flag = zeros(n,m);
%     % if dtmod2(j,i) == dtmod(j,i)
%     %    flag(j,i) = 1;
%     % end
%     % flag
%     %%%%%%%%%%%%-----------------------------
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
UU=Qnew;
CC=Cnew;
% pause
