function out=rpnew2(u,v,U,miuo,miuw,alpha,para,dt,KK)

dx=para.box.dx;
dy=para.box.dy;
m=para.box.m+1; %%%% no of columns i.e x direction of reservoir
n=para.box.n+1; %%%% no of rows i.e. y direction of reservoir

Q(:,:,1)=U(:,:,1);
Q(:,:,2)=U(:,:,1).*U(:,:,2);
%Q(:,:,3)=U(:,:,1).*U(:,:,3);


alpha0 = 0.125;
m_avg = 2/3;
% n_avg = 1.227;
temp = zeros(m,n,3);
% counter = 1;

%while(norm(Q(:,:,1)-temp(:,:,1))>10^(-3) )
%     counter = counter+1;
    temp = Q;
    %%% --redundant calculations of relative permeabilities separately---
%     kra = sqrt(Q(:,:,1)).*(1-(1-Q(:,:,1).^(1/m_avg)).^m_avg).^2;
%     kro = sqrt(Q(:,:,1)).*(1-(1-Q(:,:,1).^(1/m_avg)).^m_avg).^(2*m);
%%%%%%%%%% ------------



%     lambda_a =sqrt(Q(:,:,1)).*(1-(1-Q(:,:,1).^(1/m_avg)).^m_avg).^2 ./(miuw*(1+alpha*(Q(:,:,2)./Q(:,:,1))));
%     lambda_o =sqrt(Q(:,:,1)).*(1-(1-Q(:,:,1).^(1/m_avg)).^m_avg).^(2*m_avg)./miuo;
%     
    
    %%%%Trial----- without correction for polymer and also viscosities
    lambda_a = (Q(:,:,1).^3)./(miuw*(1+alpha*(Q(:,:,2)./Q(:,:,1))));
    lambda_o = ((1-Q(:,:,1)).^3)/miuo;
    %%%%%%% -------
    
    lambda = lambda_a+lambda_o;
    
    
    %%%%------Liqun's code-----------
%     lambda = (Q(:,:,1).^2./(miuw*(1+alpha*Q(:,:,2)./Q(:,:,1)))+(1-Q(:,:,1)).^2/miuo);
%     lambda_a = Q(:,:,1).^2./(miuw*(1+alpha*Q(:,:,2)./Q(:,:,1)));
%     lambda_o = (1-Q(:,:,1)).^2/miuo;
    %%%% parameters of capillary pressure model van Genuchten; S_w,r =0.1
    %Se = (U(:,:,1)-0.1)/(1-0.1);
    
    %%%%%%------------
    
    %%% computing derivative of P_c w.r.t s(saturation)
    Pc_diff = (1/alpha0)*((m_avg-1)/m_avg)*((Q(:,:,1).^(-1/m_avg)-1).^(-m_avg))./(Q(:,:,1).^(1+1/m_avg));
%     pause;
    %%%p_c_diff = - (m_avg*n_avg*alpha0*0.9 .*Se^((1+2*m_avg)/m_avg).*(1-Se^(1/m_avg))^(m_avg))^-1;
    %%%
    %%% computing non-linear coefficients of saturation terms
%     Dw = KK.*(lambda_a.*lambda_o./lambda).*Pc_diff;
    
    
    
    %%%%-- Trial-------
    Dw_t = 0.04*0.01.*Q(:,:,1).*(1-Q(:,:,1));
    
    
    
    %%Dc = Dw.*(Q(:,:,2)./Q(:,:,1));
  %  Dg = Dw.*(Q(:,:,3)./Q(:,:,1));
    %%% adding a pseudo computation layer of zeros around the coefficient
    %%% matrix
    %Dw_t = [zeros(1,n+2);zeros(m,1),Dw,zeros(m,1);zeros(1,n+2)]; 
    %Dc_t = [zeros(1,n+2);zeros(m,1),Dc,zeros(m,1);zeros(1,n+2)]; 
    %%%%
    %%% Computing matrix of known terms on the right hand side of the saturation equations
    Fw_ut = (lambda_a./lambda).*u;
    Fc_ut = Fw_ut.*(Q(:,:,2)./Q(:,:,1));
    Fw_vt = (lambda_a./lambda).*v;
    Fc_vt = Fw_vt.*(Q(:,:,2)./Q(:,:,1));
%     Fw_ut = [zeros(1,n+2);zeros(m,1),Fw_u,zeros(m,1);zeros(1,n+2)]; 
%     Fw_vt = [zeros(1,n+2);zeros(m,1),Fw_v,zeros(m,1);zeros(1,n+2)]; 
%     Fc_ut = [zeros(1,n+2);zeros(m,1),Fc_u,zeros(m,1);zeros(1,n+2)];
%     Fc_vt = [zeros(1,n+2);zeros(m,1),Fc_v,zeros(m,1);zeros(1,n+2)];
    %%%
    iter = 1;
    DDD = zeros((n-2)*(m-2),3);
%     AAA=zeros(m*n);
    while(iter<(n-2)*(m))
        BB=zeros(n-2,m-2);AA=zeros(n-2,m-2);CC=zeros(n-2,m-2);DD=zeros(m-2,3);
        cnt=fix(iter/(m));
        for i=2:1:n-1  %% i = x direction discretization and rows ie y direction in algorithm
            for j=2:1:m-1 %% j = y direction discretization and columns ie x direction in algorithm
                if(j==i)
                    if(iter>m)
                        AA(i-1,j-1)   = - (Dw_t(cnt+1,j)+Dw_t(cnt+2,j))/(2*dy*dy);
                        %%AA(i,j,2)   = - (Dw_t(cnt+1,j+1)+Dw_t(cnt+2,j+1))/(2*dy*dy);
                    end
                    BB(i-1,j-1)   = (.35/dt + 1/(2*dx*dx)*(Dw_t(cnt+2,j-1)+2*Dw_t(cnt+2,j)+Dw_t(cnt+2,j+1))+1/(2*dy*dy)*(Dw_t(cnt+1,j)+2*Dw_t(cnt+2,j)+Dw_t(cnt+3,j)));
                    if(iter<(n-3)*m)
                        CC(i-1,j-1)   = - (Dw_t(cnt+2,j)+Dw_t(cnt+3,j))/(2*dy*dy);
                    end
                    if(j<m-1)
                    BB(i-1,j)   = - (Dw_t(cnt+2,j)+Dw_t(cnt+2,j+1))/(2*dx*dx);
                    end
                    if(j>2)
                    BB(i-1,j-2)   = - (Dw_t(cnt+2,j)+Dw_t(cnt+2,j-1))/(2*dx*dx);
                    end
                    DD(j-1,1)     = (.35/dt)*Q(cnt+2,j,1)- (1/(2*dx))*(Fw_ut(cnt+2,j+1)-Fw_ut(cnt+2,j-1))-(1/(2*dy))*(Fw_vt(cnt+3,j)-Fw_vt(cnt+1,j));
                    DD(j-1,2)     = (.35/dt)*Q(cnt+2,j,2)- (1/(2*dx))*(Fc_ut(cnt+2,j+1)-Fc_ut(cnt+2,j-1))-(1/(2*dy))*(Fc_vt(cnt+3,j)-Fc_vt(cnt+1,j));
                end
            end
            
        end
        if(cnt==0)
           AAA = [BB,CC,zeros(m-2,(n-2)*(m-2)-2*(m-2))];
        elseif(cnt==(m-3))
           AAA = [AAA;zeros(m-2,(n-2)*(m-2)-2*(m-2)),AA,BB];
        else
           AAA = [AAA;zeros(m-2,(cnt-1)*(m-2)),AA,BB,CC,zeros(m-2,(n-2)*(m-2)-3*(m-2)-(cnt-1)*(m-2))];
        end
        
        
        DDD(cnt*(m-2)+1:cnt*(m-2)+m-2,1) = DD(:,1);
        
        DDD(cnt*(m-2)+1:cnt*(m-2)+m-2,2) = DD(:,2);
            
        iter =iter+m;
    end
%     AAA
%     size(AAA)
%     DDD
%     size(DDD)
%     Q(:,:,1)
%     pause
    
    
%     M=reshape(Q(:,:,1)',m*n,1);
    intd=(reshape(bicgstab(AAA,DDD(:,1),1*10^(-2),300,[],[]),n-2,m-2))';
    Q(:,:,1)=[Q(1,:,1);Q(2:n-1,1,1),intd,Q(2:n-1,m,1);Q(n,:,1)];
    disp('rpnew1: I am here')
    intc=(reshape(bicgstab(AAA,DDD(:,2),1*10^(-2),300),n-2,m-2))';
    disp('rpnew1: Now I am here')
    Q(:,:,2)=[Q(1,:,2);Q(2:n-1,1,2),intc,Q(2:n-1,m,2);Q(n,:,2)];
%     pause;
%end


Q(:,:,1)=Q(:,:,1);
%Q(1,1,1)=1;
Q(:,:,2)=Q(:,:,2)./Q(:,:,1);
%Q(1,1,2)=0.05;
%Q(:,:,3)=Q(:,:,3)./Q(:,:,1);

%   Q(:,1,1)=Q(:,2,1);
%   Q(1,:,1)=Q(2,:,1);
% %  
%    Q(:,end,1)=Q(:,end-1,1);
%    Q(end,:,1)=Q(end-1,:,1);
%   Q(:,end,2)=-Q(:,end,2);
%   Q(end,:,2)=-Q(end,:,2);

out=Q;