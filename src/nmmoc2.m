function UU = nmmoc2(u,v,~,~,S,C,miuo,miuw,para,dt)
% code to compute solution of saturation equations by Modified Method of
% Characteristics using explicit formulation (Yuan Yi-Rang 1993)
format long
m=size(S,2); n=size(S,1); % Notice m = para.box.m+1 and n = para.box.n+1
dx = para.box.dx; dy = para.box.dy;
Q = S;

%make an array of dt
dtcal = dt*ones(n,m);

% fixed coordinates
[x,y] = meshgrid(para.box.left:para.box.dx:para.box.right,para.box.bottom:para.box.dy:para.box.top);

% defining coefficient of hyperbolic part, b(s)
%%%% b(s) = df/ds(s) for flow without polymer
b_s = eval_bs(Q,C,miuw,miuo,0);


%determine the correct dt so that backward differencing along
%characteristic tangent does not fall outside domain

% redefined coordinates
[xmod,ymod,dtmod] = eval_X(x,y,b_s,u,v,dtcal);
xmod(:,1:5);
ymod(:,1:5);
% pause

% bilinear interpolant for saturation on redefined coordinates
Qmod = interp2(x,y,Q,xmod,ymod);
% Qmod(:,1:5);
% b_s = eval_bs(Qmod,miuw,miuo);
% % optimally redefined coordinates
% [xmod,ymod] = eval_X(x,y,b_s,umod,vmod,dt);
% % saturation values on redefined coordinates interpolated from knot values
% Qmod = interp2(x,y,Q,xmod,ymod,'linear',0);


D = 0.04*0.01.*Qmod.*(1-Qmod);
Qnew=zeros(n-2,m-2);  % initialize Qnew for speed of computation

%%%% Solve saturation equation -----------
for i=2:1:m-1  %% i = x direction discretization formula and columns ie 2nd coordinate in matrix code
    for j=2:1:n-1 %% j = y direction discretization formula and rows ie 1st coordinate in matrix code        

%%%%  Explicit formulation without polymer

       Qnew(j-1,i-1) = Qmod(j,i)+ dtmod(j,i)*(((1/(2*dx*dx))*(D(j,i-1)+2*D(j,i)+D(j,i+1))...
           +(1/(2*dy*dy))*(D(j-1,i)+2*D(j,i)+D(j+1,i)))*Q(j,i)...
           +((D(j,i)+D(j,i-1))/(2*dx*dx))*Q(j,i-1)+((D(j,i)+D(j-1,i))/(2*dy*dy))...
           *Q(j-1,i)+((D(j,i)+D(j,i+1))/(2*dx*dx))*Q(j,i+1)+((D(j,i)+D(j+1,i))...
           /(2*dy*dy))*Q(j+1,i));        

%        Qnew(j-1,i-1) = (1/(1/dtmod(j,i) + (1/(2*dx*dx))*(D(j,i-1)+2*D(j,i)+D(j,i+1))...
%            +(1/(2*dy*dy))*(D(j-1,i)+2*D(j,i)+D(j+1,i))))*((1/dtmod(j,i))*Qmod(j,i)...
%            +((D(j,i)+D(j,i-1))/(2*dx*dx))*Q(j,i-1)+((D(j,i)+D(j-1,i))/(2*dy*dy))...
%            *Q(j-1,i)+((D(j,i)+D(j,i+1))/(2*dx*dx))*Q(j,i+1)+((D(j,i)+D(j+1,i))...
%            /(2*dy*dy))*Q(j+1,i));        
    end
end
a = Qnew(1,m-2); %right bottom vertex of Q
b = Qnew(n-2,1); %left top vertex of Q
c = Qnew(n-2,m-2); %right top vertex of Q
Q = [0.79 Qnew(1,:) a; 
    Qnew(:,1) Qnew Qnew(:,m-2);
    b Qnew(n-2,:) c];
Qnew=Q;



disp('mmoc calculation done');



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

