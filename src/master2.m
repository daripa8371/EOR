%%% Main program with automatic runs on different initial polymer 
% concentrations C00 and a grid size 16X16
clc
clear all;
clear global all;
format long

% [obs,obs_header,obs_title]=read_eas('syn_perm.out');

% matlab structure array for defining mesh and grid points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
para.box.left=0;
para.box.right=1;
para.box.bottom=0;
para.box.top=1;
para.box.m=15;     % x grid points is 0 to para.box.m
para.box.n=15;     % y grid points is 0 to para.box.n
para.box.dx = (para.box.right-para.box.left)/para.box.m;
para.box.dy = (para.box.top-para.box.bottom)/para.box.n;
[x,y] = meshgrid(para.box.left:para.box.dx:para.box.right,para.box.bottom:para.box.dy:para.box.top);

phi_test=get_phi_test(para);

%Defining right hand side of elliptic system - source terms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f=zeros(para.box.n+1,para.box.m+1);
src = 200;    % magnitude of mass flow rate at source 
f(1,1)=src;                           % magnitude of injection well = 500
f(para.box.n+1,para.box.m+1)=-src;   % magnitude of production well = 500

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('KK16.mat');
figure(3)
surf(x,y,KK);
export_fig(sprintf('KK%dx%d.pdf',para.box.m+1, para.box.n+1));
%%%%%%%%%%% Preparing automatic runs

C00 = [0.1 0.01 0.001];
time = zeros(1,3);
for counter = 1:3      
    tic
    
    
    % initialize RP variables.............................................
    % UU=zeros(para.box.n+1,para.box.m+1);
    s0=0.79; c0=C00(counter);
    [UU,CC]=s0c0(para,phi_test,s0,c0);
    
    % Defining a random matrix for heterogeneous coefficient
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
%     KK=randn(para.box.n+1, para.box.m+1);
%     KK=exp(KK);
%     % KK =ones(para.box.n+1, para.box.m+1);  %%% For a constant permeability
%     
%     surf(x,y,KK);
%     export_fig(sprintf('KK%dx%dc%g.pdf',para.box.m+1, para.box.n+1, c0));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%
    miuo=12.6;  %%% 0.95
    miuw=0.0430;   %%% 0.095
    alpha=0.6;  %%% 0.6
    
    tstop=400;
    t=0;
    cfl=0.2;
    u_old = zeros(para.box.n+1,para.box.m+1); v_old = u_old;
    dt=0.01;
    
    while (t<tstop && UU(para.box.n+1,para.box.m+1) <= 0.41)
        
        
%         if t < 4
%             dt = dt0;
%         else
%             dt = 10*dt0;
%         end
        
        t=t+dt;
        
        %Solve elliptic equation.................................
        %%%%%%%% Mobilities Lambda = (k_relative/viscosity)
        %%%%Trial-----
        lambda_a = UU.^3./(miuo*(0.5+CC));
        %     lambda_a = Q(:,:,1).^3./(miuw*(1+alpha*(Q(:,:,2)./Q(:,:,1))));
        lambda_o = (1-UU).^3./miuo;
        %%%%%%% -------
        
        lambda = lambda_a +lambda_o;
        beta=KK.*lambda;  %%% for polymer flow
        %     beta=KK.*((UU.^3)./(miuw)+((1-UU).^3)./miuo);  %%% for non-polymer
        
        [U,L] = setTri(para);
        grid = setGrid(para,U,L,beta);
        rh = setRightHand(para,f,U,L);
        A = setA(para,grid);
        A = spconvert(A);
        B = setB(para,grid,rh);
        u = getu(A,B);
        vn=get_vn(u,para);
        
        %Renew Speed, phi, s, c....................................
        
        [px,py]=get_gra(vn,para);
        %     [px,py]=gradient(vn,para.box.dx,para.box.dy);
        u=-beta.*px;
        v=-beta.*py;
        %     figure(5)
        %     quiver(x,y,u,v);
        
        
        %Solve riemann problem....................................
        
        %%%%  For non-polymer flow
        %     UU=nmmoc(u,v,u_old,v_old,UU,CC,miuo,miuw,para,dt,src,s0,c0);
        %%%%  For polymer flow
        %     [UU,CC]=nmmoc3(u,v,u_old,v_old,UU,CC,miuo,miuw,para,dt);
        [UU,CC]=nmmoc5(u,v,UU,CC,miuo,miuw,para,dt,src,s0,c0);
        %     if t == dt
        %         UU1=UU; CC1= CC;
        %         [UU,CC]=nmmoc6(u,v,UU1,CC1,UU1,CC1,miuo,miuw,para,dt,src,s0,c0);
        %     else
        %         UU0 = UU1; CC0 = CC1;
        %         UU1 = UU;  CC1 = CC;
        %         [UU,CC]=nmmoc6(u,v,UU1,CC1,UU0,CC0,miuo,miuw,para,dt,src,s0,c0);
        %     end
        u_old = u;
        v_old = v;
        
        figure(1)
        contourf(x,y,UU);
        caxis manual
        caxis([0.15 0.88]);
        view(0,90);
        colorbar
        shading flat;
        
        
        tcal = uint16(t*100)
        
        if rem(tcal,50)==0
%             tcalc=tcal/10;
            export_fig(sprintf('./Run5/ST%dc%g.pdf', tcal, c0));
            %         savefig(['ST' num2str(tcal) 'c05'],'pdf');
        end
        
        %
        caxis_max= c0 + 0.5*c0;
        figure(2)
        contourf(x,y,CC);
        caxis manual
        caxis([0 caxis_max]);
        view(0,90);
        colorbar
        shading flat;
        
        if rem(tcal,50)==0
            export_fig(sprintf('./Run5/CT%dc%g.pdf', tcal, c0));
            %         savefig(['CT' num2str(tcal) 'c05'],'pdf');
        end
        %     M1(:,round(t*1000)) = getframe;
        %          pause
    end
    time(counter) =toc;
end
% movie(M1);
% movie2avi(M1,'Entry_final.avi');

