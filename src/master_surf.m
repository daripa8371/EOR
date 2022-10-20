%%% Main program with automatic runs on different initial polymer 
% concentrations C00 and a grid size 16x16
clc
clear all;
clear global all;
format long

% declaring global variables to reduce computational time in passing
% them through function calls
global miuo miuw swr0 sor0 dt KK;

% [obs,obs_header,obs_title]=read_eas('syn_perm.out');

% matlab structure array for defining mesh and grid points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
para.box.left=0;
para.box.right=1;
para.box.bottom=0;
para.box.top=1;
para.box.m=23;     % x grid points is 0 to para.box.m
para.box.n=23;     % y grid points is 0 to para.box.n
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
% Defining a random matrix for heterogeneous coefficient
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

load('KK24.mat');
%figure(3)
%surf(x,y,KK);
% export_fig(sprintf('KK%dx%d.pdf',para.box.m+1, para.box.n+1));
%%%%%%%%%%% Preparing automatic runs

% ------ effect of polymer
C00 = [0.3 0.1 0.01 0.001];

% % ------ effect of surfactant
% G00 = [0.001 0.005 0.01 0.1];
time = zeros(1,4); iterations = zeros(1,4);
for counter = 4:4      
    tic
    
    % initialize RP variables.............................................
    % UU=zeros(para.box.n+1,para.box.m+1);
    s0=0.79;  % initial residual saturation inside the reservoir = 1-s0
    c0=C00(counter);
    g0=0.005;
    [UU,CC,GG]=s0c0(para,phi_test,s0,c0,g0,1);
    
    
    %%%%%
    miuo=12.6;  %%% 0.95
    miuw=0.0430;   %%% 0.095
    alpha=0.6;  %%% 0.6
    % define initial residual saturations before critical capillary number 
    swr0 = 0.1;
    sor0 = 0.3;
    
    
    tstop=400;
    t=0;
    %cfl=0.2;
    u_old = zeros(para.box.n+1,para.box.m+1); v_old = u_old; u=u_old;v=v_old;
    dt=0.005;
    
    % setting up for movie 
    mvicntr = 1;
  
    while (t<tstop && UU(para.box.n+1,para.box.m+1) <= 0.40)
        
             
        t=t+dt;
        
        %Solve elliptic equation.................................
        
        % IFT as a function of surf concentration
        % $$ \sigma = \frac{10.001}{\Gamma +1} - 0.001 $$
        sigma = 10.001./(GG+1)-0.001;
        
        % aq soln viscosity as a function of polymer $$ \mu_a =
        % \mu_o(0.5+c) $$
        miua = compvis(CC);
        
        % recompute residual saturations $$ s_{ra}, s_{ro} $$
        [swr,sor] = compres(sigma,u,v,miua);
 
        % recompute mobilities ( with surfactants )
        lambda_a = compmob(UU,CC,sor,swr,1,1);
        lambda_o = compmob(UU,CC,sor,swr,0,1);
        
        lambda = lambda_a + lambda_o;
        
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
        %     Plot of velocity field
        %     figure(5)
        %     quiver(x,y,u,v);
        
        
        %Solve riemann problem....................................
        
        %%%%  For non-polymer flow
        %     UU=nmmoc(u,v,u_old,v_old,UU,CC,miuo,miuw,para,dt,src,s0,c0);
        %%%%  For polymer flow
        %     [UU,CC]=nmmoc3(u,v,u_old,v_old,UU,CC,miuo,miuw,para,dt);
        [UU,CC,GG]=nmmoc_surf_mod(u,v,UU,CC,GG,miua,para,sigma,c0,g0,src);

        u_old = u;
        v_old = v;
        
        tcal = uint16(t*1000); msg = ['tcal = ',num2str(tcal)]; disp(msg);
        
        saxis_max = s0 + 0.12*s0;
        caxis_max = c0 + 0.25*c0;
        gaxis_max = g0 + 0.25*g0;
        
        if rem(tcal,100)==0
            tcalc = tcal/10;
            
            figure(1)
            contourf(x,y,UU);
            caxis manual
            caxis([0.21 saxis_max]);
            view(0,90);
            colorbar
            shading flat;
            
            export_fig(sprintf('./Run18/ST%dc%g.pdf', tcalc, c0));
            M(mvicntr) = getframe(gcf);            
            
            
            figure(2)
            contourf(x,y,CC);
            caxis manual
            caxis([0 caxis_max]);
            view(0,90);
            colorbar
            shading flat;
            
            export_fig(sprintf('./Run18/CT%dc%g.pdf', tcalc, c0));
            N(mvicntr) = getframe(gcf);
            
            
            figure(3)
            contourf(x,y,GG);
            caxis manual
            caxis([0 gaxis_max]);
            view(0,90);
            colorbar
            shading flat;
            
            export_fig(sprintf('./Run18/GT%dc%g.pdf', tcalc, c0));
            P(mvicntr) = getframe(gcf);
            mvicntr = mvicntr+1;
        end
        
    end
    time(counter) =toc;
    iterations(counter) = tcal;
end
save runtime.mat time iterations
movie(M);
movie2avi(M,'Final_Sat_24_c0001.avi');
movie2avi(N,'Final_Pol_24_c0001.avi');
movie2avi(P,'Final_Surf_24_c0001.avi');


        