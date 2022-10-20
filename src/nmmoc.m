function out = nmmoc(u,v,u_old,v_old,S,~,miuo,miuw,para,dt,src,s0,c0)
% code to compute solution of saturation equations by Modified Method of
% Characteristics using implicit formulation(Yuan Yi-Rang 1993)
format long
m=size(S,2); n=size(S,1); % Notice m = para.box.m+1 and n = para.box.n+1
dx = para.box.dx; dy = para.box.dy;
Q = S;
% extrapolated velocity to (n+1)th time
umod = 2*u-u_old;
vmod = 2*v-v_old;
%make an array of dt
dtcal = dt*ones(n,m);

% fixed coordinates
[x,y] = meshgrid(para.box.left:para.box.dx:para.box.right,para.box.bottom:para.box.dy:para.box.top);

% defining coefficient b(s)
b_s = eval_bs(Q,Q,miuw,miuo,0);


%determine the correct dt so that backward differencing along
%characteristic tangent does not fall outside domain

% redefined coordinates
[xmod,ymod,dtmod] = eval_X(x,y,b_s,u,v,dtcal);

% pause

% bilinear interpolant for saturation on redefined coordinates
Qmod = interp2(x,y,Q,xmod,ymod);

b_s = eval_bs(Qmod,Qmod,miuw,miuo,0);
% optimally redefined coordinates
[xmod,ymod,dtmod] = eval_X(x,y,b_s,umod,vmod,dtmod);
% saturation values on redefined coordinates interpolated from knot values
Qmod = interp2(x,y,Q,xmod,ymod);

%%%%%%%% Mobilities Lambda = (k_relative/viscosity)
    %%%%Trial-----
    lambda_a = Q.^3./(miuw);
%     lambda_a = Q(:,:,1).^3./(miuw*(1+alpha*(Q(:,:,2)./Q(:,:,1))));
    lambda_o = (1-Q).^3./miuo;
    %%%%%%% -------
    
    lambda = lambda_a +lambda_o;

D = 0.04*0.01.*Q.*(1-Q);
%%%%debugging for vertical line problem
% D = ones(size(Qmod));
iter =1;
% AAA = zeros((n-2)*(m-2));
% DDD=zeros((n-2)*(m-2),1);
DDD = [];

while(iter<=(m)*(n-1)+1)
    cnt = (iter -1)/m; %% cnt = 0,1,2,3.... for iter = 1, m+1, 2m+1, 3m+1 etc...
    BB=zeros(n,m);AA=BB;CC=BB;DD=zeros(m,1);
    for i=1:1:m  %% i = x direction discretization formula and columns ie 2nd coordinate in matrix code
        for j=1:1:n %% j = y direction discretization formula and rows ie 1st coordinate in matrix code
            if j==i
                if iter == 1
                    if i == 1
                        DD(i) = (1/dt)*Qmod(cnt+1,i);%+(1-lambda_a(cnt+1,i)/lambda(cnt+1,i))*src;%+((D(j,i)+D(j,i-1))/(2*dx*dx))*Q(j,i-1)+((D(j,i)+D(j-1,i))/(2*dy*dy))*Q(j-1,i);
                        BB(j,i)   =  1/dt+((1/(dx*dx))*(D(cnt+1,i)-D(cnt+1,i+1))+(1/(dy*dy))*(D(cnt+1,i)-D(cnt+2,i)));
                        BB(j,i+1) =  (D(cnt+1,i+1)-D(cnt+1,i))/(dx*dx);
                        CC(j,i)   =  (D(cnt+2,i)-D(cnt+1,i))/(dy*dy);
                    elseif i == m
                        DD(i) = (1/dt)*Qmod(cnt+1,i);%+((D(j,i)+D(j,i+1))/(2*dx*dx))*Q(j,i+1)+((D(j,i)+D(j-1,i))/(2*dy*dy))*Q(j-1,i);
                        BB(j,i-1) =  (D(cnt+1,i-1)-D(cnt+1,i))/(dx*dx);
                        BB(j,i)   =  1/dt+((1/(dx*dx))*(D(cnt+1,i)-D(cnt+1,i-1))+(1/(dy*dy))*(D(cnt+1,i)-D(cnt+2,i)));
                        CC(j,i)   =  (D(cnt+2,i)-D(cnt+1,i))/(dy*dy);
                    else
                        DD(i) = (1/dt)*Qmod(cnt+1,i);%+((D(j,i)+D(j-1,i))/(2*dy*dy))*Q(j-1,i);
                        BB(j,i)   = 1/dt-((1/(2*dx*dx))*(D(cnt+1,i-1)+D(cnt+1,i+1))+(1/(dy*dy))*(D(cnt+2,i)));
                        BB(j,i-1) =  (D(cnt+1,i-1)+D(cnt+1,i))/(2*dx*dx);
                        BB(j,i+1) =  (D(cnt+1,i)+D(cnt+1,i+1))/(2*dx*dx);
                        CC(j,i)   =  (D(cnt+2,i)-D(cnt+1,i))/(dy*dy);
                    end
                elseif iter == (m)*(n-1)+1
                    if i == 1
                        DD(i) = (1/dt)*Qmod(cnt+1,i);%+((D(j,i)+D(j,i-1))/(2*dx*dx))*Q(j,i-1)+((D(j,i)+D(j+1,i))/(2*dy*dy))*Q(j+1,i);
                        AA(j,i)   =  (D(cnt,i)-D(cnt+1,i))/(dy*dy);
                        BB(j,i)   =  1/dt+((1/(dx*dx))*(D(cnt+1,i)-D(cnt+1,i+1))+(1/(dy*dy))*(D(cnt+1,i)-D(cnt,i)));
                        BB(j,i+1) =  (D(cnt+1,i+1)-D(cnt+1,i))/(dx*dx);
                    elseif i == m
                        DD(i) = (1/dt)*Qmod(cnt+1,i);%+((D(j,i)+D(j,i+1))/(2*dx*dx))*Q(j,i+1)+((D(j,i)+D(j+1,i))/(2*dy*dy))*Q(j+1,i);
                        AA(j,i)   =  (D(cnt,i)-D(cnt+1,i))/(dy*dy);
                        BB(j,i-1) =  (D(cnt+1,i-1)-D(cnt+1,i))/(dx*dx);
                        BB(j,i)   = 1/dt-((1/(dx*dx))*(D(cnt+1,i-1)-D(cnt+1,i))+(1/(dy*dy))*(D(cnt,i)-D(cnt+1,i)));
                    else
                        DD(i) = (1/dt)*Qmod(cnt+1,i);%+((D(j,i)+D(j+1,i))/(2*dy*dy))*Q(j+1,i);
                        AA(j,i)   =  (D(cnt,i)-D(cnt+1,i))/(dy*dy);
                        BB(j,i)   = 1/dt-((1/(2*dx*dx))*(D(cnt+1,i-1)+D(cnt+1,i+1))+(1/(dy*dy))*(D(cnt,i)));
                        BB(j,i-1) =  (D(cnt+1,i-1)+D(cnt+1,i))/(2*dx*dx);
                        BB(j,i+1) =  (D(cnt+1,i)+D(cnt+1,i+1))/(2*dx*dx);
                    end
                else
                    if i == 1
                        DD(i) = (1/dt)*Qmod(cnt+1,i);%+((D(j,i)+D(j,i-1))/(2*dx*dx))*Q(j,i-1);
                        AA(j,i)   =  (D(cnt,i)+D(cnt+1,i))/(2*dy*dy);
                        BB(j,i)   = 1/dt-((1/(dx*dx))*(D(cnt+1,i+1))+(1/(2*dy*dy))*(D(cnt,i)+D(cnt+2,i)));
                        BB(j,i+1) =  (D(cnt+1,i+1)-D(cnt+1,i))/(dx*dx);
                        CC(j,i)   =  (D(cnt+1,i)+D(cnt+2,i))/(2*dy*dy);
                    elseif i == m
                        DD(i) = (1/dt)*Qmod(cnt+1,i);%+((D(j,i)+D(j,i+1))/(2*dx*dx))*Q(j,i+1);
                        AA(j,i)   =  (D(cnt,i)+D(cnt+1,i))/(2*dy*dy);
                        BB(j,i-1) =  (D(cnt+1,i-1)-D(cnt+1,i))/(dx*dx);
                        BB(j,i)   = 1/dt-((1/(dx*dx))*(D(cnt+1,i-1))+(1/(2*dy*dy))*(D(cnt,i)+D(cnt+2,i)));
                        CC(j,i)   =  (D(cnt+1,i)+D(cnt+2,i))/(2*dy*dy);
                    else
                        DD(i) = (1/dt)*Qmod(cnt+1,i);
                        AA(j,i)   =  (D(cnt+1-1,i)-D(cnt+1,i))/(dy*dy);
                        BB(j,i-1) =  (D(cnt+1,i-1)-D(cnt+1,i))/(dx*dx);
                        BB(j,i)   = 1/dt+((1/(dx*dx))*(D(cnt+1,i)-D(cnt+1,i-1))+...
                            (1/(dy*dy))*(D(cnt+1,i)-D(cnt,i)));
%                         AA(j,i)   =  (D(cnt+1-1,i)+D(cnt+1,i))/(2*dy*dy);
%                         CC(j,i)   =  (D(cnt+1,i)+D(cnt+1+1,i))/(2*dy*dy);
%                         BB(j,i)   = 1/dt-((1/(2*dx*dx))*(D(cnt+1,i-1)+2*D(cnt+1,i)+D(cnt+1,i+1))+...
%                             (1/(2*dy*dy))*(D(cnt+1-1,i)+2*D(cnt+1,i)+D(cnt+1+1,i)));
%                         BB(j,i+1) =  (D(cnt+1,i)+D(cnt+1,i+1))/(2*dx*dx);
%                         BB(j,i-1) =  (D(cnt+1,i-1)+D(cnt+1,i))/(2*dx*dx);
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
% Q = reshape(QQQ,m,n)';

Qnew=(reshape(bicgstab(AAA,DDD,1*10^(-10),600,[],[]),m,n))';
% Qnew(1,1)=s0;
disp('mmoc calculation done');
% Q(:,1:5);

% a = Q(1,m-2); %right bottom vertex of Q
% b = Q(n-2,1); %left top vertex of  Q
% c = Q(n-2,m-2); %right top vertex of Q
% Q = [0.79 Q(1,:) a; 
%     Q(:,1) Q Q(:,m-2);
%     b Q(n-2,:) c];
Qnew(:,1:5)
% pause;
% Qnew=S;
% Qnew(2:n-1,2:m-1) = Q;

figure(100)
plot(x(1,:),S(2,:),'Color','blue');
hold on
axis([0 1 -0.5 0.9])
plot(x(1,:),Qnew(2,:)-0.4,'Color','red');


% figure(101)
% plot(x(1,:),S(floor(n-2),:),'Color','blue');
% hold on
% axis([0 1 -0.5 0.9])
% plot(x(1,:),Qnew(floor(n-2),:)-0.4,'Color','red');
out=Qnew;

