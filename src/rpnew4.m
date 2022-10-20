function out=rpnew4(u,v,U,miuo,miuw,alpha,para,dt,KK,src)

Q(:,:,1)=U(:,:,1);
%Q(:,:,2)=U(:,:,1).*U(:,:,2);
%Q(:,:,3)=U(:,:,1).*U(:,:,3);
dx=para.box.dx;
dy=para.box.dy;
m=para.box.m + 1;
n=para.box.n + 1;


%alpha0 = 0.125;
%m_avg = 2/3;
% n_avg = 1.227;
temp = zeros(size(Q,1),size(Q,2),3);
% counter = 1;
error=[];


while(abs(norm(Q(:,:,1)-temp(:,:,1)))>10^(-4))
     error = [error,abs(norm(Q(:,:,1)-temp(:,:,1)))];
     temp = Q;

    %%% --redundant calculations of relative permeabilities separately---
%     kra = sqrt(Q(:,:,1)).*(1-(1-Q(:,:,1).^(1/m_avg)).^m_avg).^2;
%     kro = sqrt(Q(:,:,1)).*(1-(1-Q(:,:,1).^(1/m_avg)).^m_avg).^(2*m);
%%%%%%%%%% ------------



% %     lambda_a =sqrt(Q(:,:,1)).*(1-(1-Q(:,:,1).^(1/m_avg)).^m_avg).^2 ./(miuw*(1+alpha*(Q(:,:,2)./Q(:,:,1))));
% %     lambda_o =sqrt(Q(:,:,1)).*(1-(1-Q(:,:,1).^(1/m_avg)).^m_avg).^(2*m_avg)./miuo;
    
    %%%%%%%% Mobilities Lambda = (k_relative/viscosity)
    %%%%Trial-----
    lambda_a = temp(:,:,1).^3./(miuw);
%     lambda_a = Q(:,:,1).^3./(miuw*(1+alpha*(Q(:,:,2)./Q(:,:,1))));
    lambda_o = (1-temp(:,:,1)).^3./miuo;
    %%%%%%% -------
    
    lambda = lambda_a +lambda_o;
    
    
    %%%%------Liqun's code-----------
%     lambda = (U(:,:,1).^2./(miuw*(1+alpha*U(:,:,2)))+(1-U(:,:,1)).^2/miuo);
%     lambda_a = U(:,:,1).^2./(miuw*(1+alpha*U(:,:,2)));
%     lambda_o = (1-U(:,:,1)).^2/miuo;
    %%%% parameters of capillary pressure model van Genuchten; S_w,r =0.1
    %Se = (U(:,:,1)-0.1)/(1-0.1);
    
    %%%%%%------------
    
    %%% computing derivative of P_c w.r.t s(saturation)
%    Pc_diff = (1/alpha0)*((m_avg-1)/m_avg)*((Q(:,:,1).^(-1/m_avg)-1).^(-m_avg))./(Q(:,:,1).^(1+1/m_avg));
%     pause;
    %%%p_c_diff = - (m_avg*n_avg*alpha0*0.9 .*Se^((1+2*m_avg)/m_avg).*(1-Se^(1/m_avg))^(m_avg))^-1;
    %%%
  
    %%% Computing matrix of known terms on the right hand side of the saturation equations
    Fw_u = (lambda_a./lambda).*u;
%     Fc_u = Fw_u.*(Q(:,:,2)./Q(:,:,1));
    Fw_v = (lambda_a./lambda).*v;
%     Fc_v = Fw_v.*(Q(:,:,2)./Q(:,:,1));
 %   Fw_ut = [zeros(1,size(Q,2)+2);zeros(size(Q,1),1),Fw_u,zeros(size(Q,1),1);zeros(1,size(Q,2)+2)]; 
  %  Fw_vt = [zeros(1,size(Q,2)+2);zeros(size(Q,1),1),Fw_v,zeros(size(Q,1),1);zeros(1,size(Q,2)+2)]; 
%     Fc_ut = [zeros(1,size(Q,2)+2);zeros(size(Q,1),1),Fc_u,zeros(size(Q,1),1);zeros(1,size(Q,2)+2)];
%     Fc_vt = [zeros(1,size(Q,2)+2);zeros(size(Q,1),1),Fc_v,zeros(size(Q,1),1);zeros(1,size(Q,2)+2)];
    

    for i=2:n-1
        for j=2:m-1
            if(i==2 && j ==2)
                Q(i,j,1) = U(i,j,1)-(dt/(2*dx))*(Fw_u(i,j+1)-Fw_u(i,j-1))-(dt/(2*dy))*(Fw_v(i+1,j)-Fw_v(i-1,j))+dt*(1-lambda_a(i,j)/lambda(i,j))*src;
            else
                Q(i,j,1) = U(i,j,1)-(dt/(2*dx))*(Fw_u(i,j+1)-Fw_u(i,j-1))-(dt/(2*dy))*(Fw_v(i+1,j)-Fw_v(i-1,j));

%             if(i==1 && j ==1)
%                 Q(i,j,1) = U(i,j,1)-(dt/(2*dx))*(Fw_ut(i+1,j+2)-Fw_ut(i+1,j))-(dt/(2*dy))*(Fw_vt(i+2,j+1)-Fw_vt(i,j+1))+(1-lambda_a(i,j)/lambda(i,j))*src;
%             else
%                 Q(i,j,1) = U(i,j,1)-(dt/(2*dx))*(Fw_ut(i+1,j+2)-Fw_ut(i+1,j))-(dt/(2*dy))*(Fw_vt(i+2,j+1)-Fw_vt(i,j+1));
%               DD(j,2)  = (2/dt)*Q(cnt+1,j,2)- (1/(2*dx))*(Fc_ut(cnt+2,j+2)-Fc_ut(cnt+2,j))-(1/(2*dy))*(Fc_vt(cnt+3,j+1)-Fc_vt(cnt+1,j+1));
            end
        end
    end
%     while(iter<=size(Q,1)*size(Q,2))
%         BB=zeros(size(Q,1),size(Q,2));AA=zeros(size(Q,1),size(Q,2));CC=zeros(size(Q,1),size(Q,2));DD=zeros(size(Q,2),3);
%         cnt=fix(iter/size(Q,2));
%         for i=1:1:size(Q,1)  %% i = x direction discretization formula and rows ie y direction in code
%             for j=1:1:size(Q,2) %% j = y direction discretization formula and columns ie x direction in code
%                 if(j==i)
%                     if(iter>size(Q,2))
%                         %AA(i,j,1)   = - (Dw(cnt+1,j+1)+Dw(cnt+2,j+1))/(2*dy*dy);
%                         AA(i,j)   = - (Dw_t(cnt+1,j+1)+Dw_t(cnt+2,j+1))/(2*dy*dy);
% %                         AA(i,j,2)   = - (Dw_t(cnt+1,j+1)+Dw_t(cnt+2,j+1))/(2*dy*dy);
%                     end
%                     %BB(i,j)   = (1/dt + 1/(2*dx*dx)*(Dw(cnt+2,j)+2*Dw(cnt+2,j+1)+Dw(cnt+2,j+2))+1/(2*dy*dy)*(Dw(cnt+1,j+1)+2*Dw(cnt+2,j+1)+Dw(cnt+3,j+1)));
%                     BB(i,j)   = (1/dt + 1/(2*dx*dx)*(Dw_t(cnt+2,j)+2*Dw_t(cnt+2,j+1)+Dw_t(cnt+2,j+2))+1/(2*dy*dy)*(Dw_t(cnt+1,j+1)+2*Dw_t(cnt+2,j+1)+Dw_t(cnt+3,j+1)));
%                     if(iter<size(Q,1)*size(Q,2))
%                         CC(i,j)   = - (Dw_t(cnt+2,j+1)+Dw_t(cnt+3,j+1))/(2*dy*dy);
%                     end
%                     if(j<size(Q,2))
%                     BB(i,j+1)   = - (Dw_t(cnt+2,j+1)+Dw_t(cnt+2,j+2))/(2*dx*dx);
%                     end
%                     if(j>1)
%                     BB(i,j-1)   = - (Dw_t(cnt+2,j+1)+Dw_t(cnt+2,j))/(2*dx*dx);
%                     end
%                     if(iter==0 && j==1)
%                         DD(j,1)     = (1/dt)*U(cnt+1,j,1)- (1/(2*dx))*(Fw_ut(cnt+2,j+2)-Fw_ut(cnt+2,j))-(1/(2*dy))*(Fw_vt(cnt+3,j+1)-Fw_vt(cnt+1,j+1))+(1-lambda_a(j,j)/lambda(j,j))*src;
%                     else
%                         DD(j,1)     = (1/dt)*U(cnt+1,j,1)- (1/(2*dx))*(Fw_ut(cnt+2,j+2)-Fw_ut(cnt+2,j))-(1/(2*dy))*(Fw_vt(cnt+3,j+1)-Fw_vt(cnt+1,j+1));
% %                     DD(j,2)     = (2/dt)*Q(cnt+1,j,2)- (1/(2*dx))*(Fc_ut(cnt+2,j+2)-Fc_ut(cnt+2,j))-(1/(2*dy))*(Fc_vt(cnt+3,j+1)-Fc_vt(cnt+1,j+1));
%                     end
%                 end
%             end
%             
%         end
%         if(iter==1)
%            AAA = [BB,CC,zeros(size(Q,1),size(Q,1)*size(Q,2)-2*size(Q,2))];
%         elseif(rem(iter-1,size(Q,2))==0 && (iter-1)/size(Q,2)<size(Q,1)-1)
%            AAA = [AAA;zeros(size(Q,1),(iter-1)-size(Q,2)),AA,BB,CC,zeros(size(Q,1),size(Q,1)*size(Q,2)-3*size(Q,2)-(iter-1)+size(Q,2))];
%         elseif(iter==(size(Q,1)-1)*size(Q,2)+1)
%            AAA = [AAA;zeros(size(Q,1),(iter-1+size(Q,2))-2*size(Q,2)),AA,BB];
%         end
%         
%         
%         DDD(iter:iter+size(Q,2)-1,1) = DD(:,1);
%         
% %         DDD(iter:iter+size(Q,2)-1,2) = DD(:,2);
%             
%         iter =iter+size(Q,2);
%     size(AAA)
%     DDD
%     size(DDD)
%     Q(:,:,1)
%     pause
    
    
%     M=reshape(Q(:,:,1)',size(Q,1)*size(Q,2),1);
    %Q(:,:,1)=reshape(bicgstab(AAA,DDD(:,1),1*10^(-6),600,[],[]),size(Q,1),size(Q,2))';
    
%     Q(:,:,2)=reshape(bicgstab(AAA,DDD(:,2),1*10^(-2),600,[],[]),size(Q,1),size(Q,2))';
    disp('rpnew3: I am here')
%     pause;
Q(1,1,1)=.79;
 Q(:,1,1)=Q(:,2,1);
 Q(1,:,1)=Q(2,:,1);
 Q(:,end,1)=Q(:,end-1,1);
 Q(end,:,1)=Q(end-1,:,1);
end
error


% Q(:,:,1)=Q(:,:,1);
%Q(1,1,1)=.79;
% Q(:,:,2)=Q(:,:,2)./Q(:,:,1);
% Q(1,1,2)=0.05;
%Q(:,:,3)=Q(:,:,3)./Q(:,:,1);

%    Q(:,1,1)=Q(:,2,1);
%    Q(1,:,1)=Q(2,:,1);
%    Q(:,end,1)=Q(:,end-1,1);
%    Q(end,:,1)=Q(end-1,:,1);
% Q(:,end,2)=-Q(:,end,2);
% Q(end,:,2)=-Q(end,:,2);

out=Q;