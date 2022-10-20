function out = mmoc(u,v,u_old,v_old,S,C,miuo,miuw,para,dt,~,~,~)
% code to compute solution of saturation equations by Modified Method of
% Characteristics (Douglas 1983)
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

% redefined coordinates
[xmod,ymod,dtmod] = eval_X(x,y,b_s,umod,vmod,dtcal);

% bilinear interpolant for saturation on redefined coordinates
Qmod = interp2(x,y,Q,xmod,ymod,'linear',0);

b_s = eval_bs(Qmod,Qmod,miuw,miuo,0);
% optimally redefined coordinates
[xmod,ymod,dtmod] = eval_X(x,y,b_s,umod,vmod,dtcal);
% saturation values on redefined coordinates interpolated from knot values
Qmod = interp2(x,y,Q,xmod,ymod,'linear',0);


% D = 0.04*0.01.*Qmod.*(1-Qmod);
%%%%debugging for vertical line problem
D = ones(size(Qmod)); %% setting diffusion coefficient, D(s) = 1
iter =1;
DDD=[];

while(iter<=(m)*(n-3)+1)
    cnt = (iter -1)/m; %% cnt = 0,1,2,3.... for iter = 1, m+1, 2m+1, 3m+1 etc...
    BB=zeros(n-2,m-2);AA=BB;CC=BB;DD=zeros(m-2,1);
    for i=2:1:m-1  %% i = x direction discretization formula and columns ie 2nd coordinate in matrix code
        for j=2:1:n-1 %% j = y direction discretization formula and rows ie 1st coordinate in matrix code
            if j==i
                if iter>=m+1
                    AA(j-1,i-1)   = - (D(cnt+2-1,i)+D(cnt+2,i))/(2*dy*dy);
                end
                
                if iter<(m)*(n-3)
                    CC(j-1,i-1)   = - (D(cnt+2,i)+D(cnt+2+1,i))/(2*dy*dy);
                end
                
                BB(j-1,i-1)   = 1/dt+((1/(2*dx*dx))*(D(cnt+2,i-1)+2*D(cnt+2,i)+D(cnt+2,i+1))+(1/(2*dy*dy))*(D(cnt+2-1,i)+2*D(cnt+2,i)+D(cnt+2+1,i)));
                if i<m-1
                    BB(j-1,i)   = - (D(cnt+2,i)+D(cnt+2,i+1))/(2*dx*dx);
                end
                if(i>2)
                    BB(j-1,i-2)   = - (D(cnt+2,i-1)+D(cnt+2,i))/(2*dx*dx);
                end
                if iter == 1
                    if i == 2
                        DD(i-1) = (1/dt)*Qmod(cnt+2,i)+((D(cnt+2,i)+D(cnt+2,i-1))/(2*dx*dx))*Q(j,i-1)+((D(cnt+2,i)+D(cnt+2-1,i))/(2*dy*dy))*Q(j-1,i);
                    elseif i == m-1
                        DD(i-1) = (1/dt)*Qmod(cnt+2,i)+((D(cnt+2,i)+D(cnt+2,i+1))/(2*dx*dx))*Q(j,i+1)+((D(cnt+2,i)+D(cnt+2-1,i))/(2*dy*dy))*Q(j-1,i);
                    else
                        DD(i-1) = (1/dt)*Qmod(cnt+2,i)+((D(cnt+2,i)+D(cnt+2-1,i))/(2*dy*dy))*Q(j-1,i);
                    end
                elseif iter == (m)*(n-3)+1
                    if i == 2
                        DD(i-1) = (1/dt)*Qmod(cnt+2,i)+((D(cnt+2,i)+D(cnt+2,i-1))/(2*dx*dx))*Q(j,i-1)+((D(cnt+2,i)+D(cnt+2+1,i))/(2*dy*dy))*Q(j+1,i);
                    elseif i == m-1
                        DD(i-1) = (1/dt)*Qmod(cnt+2,i)+((D(cnt+2,i)+D(cnt+2,i+1))/(2*dx*dx))*Q(j,i+1)+((D(cnt+2,i)+D(cnt+2-1,i))/(2*dy*dy))*Q(j-1,i);
                    else
                        DD(i-1) = (1/dt)*Qmod(cnt+2,i)+((D(cnt+2,i)+D(cnt+2+1,i))/(2*dy*dy))*Q(j+1,i);
                    end
                else
                    if i == 2
                        DD(i-1) = (1/dt)*Qmod(cnt+2,i)+((D(cnt+2,i)+D(cnt+2,i-1))/(2*dx*dx))*Q(j,i-1);
                    elseif i == m-1
                        DD(i-1) = (1/dt)*Qmod(cnt+2,i)+((D(cnt+2,i)+D(cnt+2,i+1))/(2*dx*dx))*Q(j,i+1);
                    else
                        DD(i-1) = (1/dt)*Qmod(cnt+2,i);
                    end
                end
                
            end
        end
    end
    if(iter==1)
        AAA = [BB,CC,zeros(n-2,(n-2)*(m-2)-2*(m-2))];
    elseif(iter==(n-3)*(m)+1)
        AAA = [AAA;zeros(n-2,(n-2)*(m-2)-2*(m-2)),AA,BB];
    else
        AAA = [AAA;zeros(n-2,iter-(m+2*cnt-1)),AA,BB,CC,zeros(n-2,(n-2)*(m-2)-3*(m-2)-(iter-(m+2*cnt-1)))];
    end
    DDD = [DDD;DD];
    iter = iter+m;
end

Q=reshape(bicgstab(AAA,DDD,1*10^(-6),600,[],[]),m-2,n-2)';
disp('mmoc calculation done');

a = Q(1,m-2); %right bottom vertex of Q
b = Q(n-2,1); %left top vertex of Q
c = Q(n-2,m-2); %right top vertex of Q
Q = [0.79 Q(1,:) a; Q(:,1) Q Q(:,m-2);b Q(n-2,:) c];

Q(:,1:5)

figure(101)
plot(x(1,:),S(floor(n-2),:),'Color','blue');
hold on
axis([0 1 -0.5 0.9])
plot(x(1,:),Q(floor(n-2),:)-0.4,'Color','red');
out=Q;

