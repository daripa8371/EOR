function out=rp(u,v,U,miuo,miuw,alpha,dx,dy,dt)

Q(:,:,1)=U(:,:,1);
Q(:,:,2)=U(:,:,1).*U(:,:,2);
%Q(:,:,3)=U(:,:,1).*U(:,:,3);

Uxp1=[U(:,2:end,:),U(:,1,:)];
Uxp2=[Uxp1(:,2:end,:),Uxp1(:,1,:)];
Uxm1=[U(:,end,:),U(:,1:end-1,:)];
Uxm2=[Uxm1(:,end,:),Uxm1(:,1:end-1,:)];

Uyp1=[U(2:end,:,:);U(1,:,:)];
Uyp2=[Uyp1(2:end,:,:);Uyp1(1,:,:)];
Uym1=[U(end,:,:);U(1:end-1,:,:)];
Uym2=[Uym1(end,:,:);Uym1(1:end-1,:,:)];

Q=Q+dt/dx*(FX(Uxm2,Uxm1,U,Uxp1,u,miuw,miuo,alpha)-FX(Uxm1,U,Uxp1,Uxp2,u,miuw,miuo,alpha))+dt/dy*(GY(Uym2,Uym1,U,Uyp1,v,miuw,miuo,alpha)-GY(Uym1,U,Uyp1,Uyp2,v,miuw,miuo,alpha));

% Q=dt*source_term(Q);

U(:,:,1)=Q(:,:,1);
U(:,:,2)=Q(:,:,2)./Q(:,:,1);
%U(:,:,3)=Q(:,:,3)./Q(:,:,1);

U(:,1,:)=U(:,2,:);
U(1,:,:)=U(2,:,:);
U(:,end,:)=U(:,end-1,:);
U(end,:,:)=U(end-1,:,:);
% U(:,end,2)=-U(:,end,2);
% U(end,:,2)=-U(end,:,2);


out=U;
