%% Main program for a full flooding simulation
% Grid Sizes used 15x15, 20x20, 30x30, 40x40 
clc
clear variables;
clear global;
format long
profile clear

%profile on
%parpool

% declaring global variables to reduce computational time in passing
% them through function calls
global miuo miuw swr0 sor0 dt KK s0 c0 g0 theta init_front_hs;

% [obs,obs_header,obs_title]=read_eas('syn_perm.out');

%%%% matlab structure array for defining mesh and grid points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
para.box.left=0;
para.box.right=1;
para.box.bottom=0;
para.box.top=1;


%%%% Initialization    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% nsim = number of runs
%%% sizeofgrid = mesh size for each run
%%% c0iter = injection concentration of polymer for each run
%%% g0iter = injection concentration of surfactant for each run
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nsim = 1;
sizeofgrid_x = [72 36 18 9]; 
sizeofgrid_y = [72 36 18 9]; 
mfine = sizeofgrid_x(1)+1; nfine=sizeofgrid_y(1)+1; % sizeofgrid is "h"
c0iter = [0.1 0 0 0];%[0.1 0.1 0.1];  %
g0iter = [0.01 0 0 0]; %[0.01 0.01 0.01];
dt = 1/50;
tstop = 100;

%%%%%
miuo=12.6;  %%% 0.95
miuw=1.26;   %%% 0.095   % Effective aq. viscosity = miuo*(0.5+c)
%alpha=0.6;  %%% 0.6
% define initial residual saturations before critical capillary number
swr0 = 0.1; % actual initial residual saturation is higher than this value in practice
sor0 = 0.2;
theta = 2;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% preallocating variables that will store cpu time taken for 
%%% each run and number of iterations/time-steps of each run
time = zeros(1,nsim); iterations = zeros(1,nsim);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%% Creating VideoWriter objects for each simulation----------
% vid1 = VideoWriter('Video1.avi'); vid1.Quality=100;
% vid2 = VideoWriter('Video2.avi'); vid2.Quality=100;
% vid3 = VideoWriter('Video3.avi'); vid3.Quality=100;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% % %%%%% Error calculation initialization -------------------------
% Smax_1 = 0; Smax_2 = 0; Smax_3 = 0; pmax_1 = 0; pmax_2 = 0; pmax_3 = 0; umax_1 = 0; umax_2 = 0; umax_3 = 0; 
% S2_c1 = 0; S2_c2 = 0; S2_c3 = 0; p2_c1 = 0; p2_c2 = 0; p2_c3 = 0; u2_c1 = 0; u2_c2 = 0; u2_c3 = 0; 
% nFrames = floor(tstop/dt);
% UU_fine = zeros(sizeofgrid_y(1)+1,sizeofgrid_x(1)+1,nFrames);
% p_fine = UU_fine; U_fine = UU_fine; V_fine = UU_fine;
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for counter = 1:nsim      
    tic
    
    para.box.m  = sizeofgrid_x(counter);     % x grid points is 0 to para.box.m
    para.box.n  = sizeofgrid_y(counter);     % y grid points is 0 to para.box.n
    para.box.dx = (para.box.right-para.box.left)/para.box.m;
    para.box.dy = (para.box.top-para.box.bottom)/para.box.n;
    [x,y] = meshgrid(para.box.left:para.box.dx:para.box.right,para.box.bottom:para.box.dy:para.box.top);

    %%% matrix containing the evaluation of the distance based level set
    % function at each of the grid points
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % WHEN FLAG = 1 it is a quarter five spot flood
    % WHEN FLAG = 2 it is a HS flood
    FLAG = 2;
    init_front_hs = 0.1;
    phi_test=get_phi_test(para,FLAG); 
        
    %Defining right hand side of elliptic system - source terms
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    f=zeros(para.box.n+1,para.box.m+1);
    if FLAG == 1
        % FOR Quarter five spot FLOOD
        src = 200;    % magnitude of mass flow rate at source
        f(1,1)=src;                           % intensity of injection well = src
        f(para.box.n+1,para.box.m+1)=-src;    % intensity of production well = -src
    elseif FLAG == 2
        
        % FOR HS FLOOD
        src = 50;    % magnitude of mass flow rate at source
        f(:,1)=src;                           % intensity of injection well = src
        f(:,para.box.m+1)=-src;    % intensity of production well = -src
        %-------------------------------------------------------------------------
    end
    
    
    %%%%%%%   Defining permeability matrix    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % flag values:  1 = homogeneous with magnitude 1000
    %               2 = continuous heterogeneous function
    %               3 = impermeable block inclusion at the center
    %               4 = impermeable block inclusion off the center
    %               5 = Upper Ness from SPE10 model sections
    %               6 = Tabert from SPE10 model sections 
    %               7 = Different Tarbert from SPE10 model sections
    %               8 = Different Ness from SPE10 model sections     
    %               9 = Heterogeneous perm function from 2016 paper 
    %               10 = Perm field with impermeable block channels
    %               11 = Permeability function from Chueh, Djilali, Bangerth 2013
    %               12 = Multiscale permeability field (32x32)
    %               13 = Grid Orientation Perm field-1
    %               14 = Grid Orientation Perm field-2
    
    KKdef(counter,para,x,y,9); 

    
%     KKdef(counter,sizeofgrid_x,sizeofgrid_y,x,y,1);
%     %%%% Visualization of Permeability field    
%     surf(x,y,KK); shading interp; colorbar; %title('Tarbert formation')
%     export_fig 'KK30x30Tarbert.pdf' '-opengl'

%     %%%% Permeability for sequential runs
%     if counter == 1
%         KKdef(counter,sizeofgrid,x,y,9)
%     elseif counter == 2
%         KKdef(counter,sizeofgrid,x,y,9)
%     end
    
    %%%%%%%%%%% Preparing automatic runs  %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % initialize RP variables----------------------------------------------
    s0=0.79;  % initial residual water saturation inside the reservoir = 1-s0
    c0=c0iter(counter);%0.1;
    g0=g0iter(counter);%0.005;
    
    [UU,CC,GG]=s0c0(para,phi_test,s0,c0,g0,logical(g0));    
    
%     fprintf('UU(1,1) = %12.10f\n',UU(1,1))
%     fprintf('UU(n,m) = %12.10f\n',UU(para.box.n+1,para.box.m+1))

%     figure(1)
%     contourf(x,y,UU,5); shading flat; view(0,90); % view for contourf plots
%     caxis manual; caxis([0.21 1.02]);colorbar; colormap jet; drawnow
    
    
    t=0; tcal = 0;
    %cfl=0.2;
    u = zeros(para.box.n+1,para.box.m+1); v = u; 
%     u_old = zeros(para.box.n+1,para.box.m+1); v_old=u_old; u=u_old; v=v_old;

    %%%% Oil recovery calculation initialization -----------------------
    COC = zeros(nsim,floor(tstop/dt));
    CROIP = zeros(nsim,floor(tstop/dt));    
    %%%% -------------------------------------
    

% %   %%%%Opening Videowriter object for writing framedata...................
%     if counter == 1
%         open(vid1);
%     elseif counter == 2
%         open(vid2);
%     elseif counter == 3
%         open(vid3);
%     end
%     %%%%%%%%%%%%%%%%%%%%%%%-------------------------------------    
    
    nstop = 3*2^(4-counter)+1;
    while (t<tstop && UU(nstop,nstop) <= 0.50)  %UU(floor(para.box.n+1/4),floor(para.box.m+1/4))
         
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
        lambda_a = compmob(UU,CC,sor,swr,1);
        lambda_o = compmob(UU,CC,sor,swr,0);        
        lambda = lambda_a + lambda_o;        
        beta = KK.*lambda;  %%% for polymer flow
        %     beta=KK.*((UU.^3)./(miuw)+((1-UU).^3)./miuo);  %%% for non-polymer
        
        [U,L] = setTri(para);
        gridsize = setGrid(para,U,L,beta);
        rh = setRightHand(para,f,U,L);
        A = setA(para,gridsize);
        A = spconvert(A);
        B = setB(para,gridsize,rh);
        p = getu(A,B); %Pressure stored as a 1d array
        vn = get_vn(p,para); %pressure stored as a nxn matrix
        
        %Renew Speed, phi, s, c....................................        
        [px,py] = get_gra(vn,para);
        %     [px,py]=gradient(vn,para.box.dx,para.box.dy);
        u = -beta.*px;
        v = -beta.*py;
        clear px py A B
        %     Plot of velocity field
        %     figure(5)
        %     quiver(x,y,u,v);
        
        %%%% Second order neumann boundary implementation
        [UU,CC,GG,OC,WC,ROIP] = nmmoc_surf_mod_neumann(u,v,UU,CC,GG,miua,para,sigma,src);

        %%%% Sequential floods  -------------------------------------------
%         if counter == 1
%             if tcal <= 250
%                 [UU,CC,GG,OC,WC,ROIP]=nmmoc_surf_mod_neumann(u,v,UU,CC,GG,miua,para,sigma,c0,g0,src);
%             else
%                 [UU,CC,GG,OC,WC,ROIP]=nmmoc_surf_mod_neumann(u,v,UU,CC,GG,miua,para,sigma,c0,0,src);
%             end
%         elseif counter ==2
%             if tcal <= 600
%                 [UU,CC,GG,OC,WC,ROIP]=nmmoc_surf_mod_neumann(u,v,UU,CC,GG,miua,para,sigma,c0,g0,src);
%             else
%                 [UU,CC,GG,OC,WC,ROIP]=nmmoc_surf_mod_neumann(u,v,UU,CC,GG,miua,para,sigma,c0,0,src);
%             end
%         elseif counter == 3
%             if tcal <= 1000
%                 [UU,CC,GG,OC,WC,ROIP]=nmmoc_surf_mod_neumann(u,v,UU,CC,GG,miua,para,sigma,c0,g0,src);
%             else
%                 [UU,CC,GG,OC,WC,ROIP]=nmmoc_surf_mod_neumann(u,v,UU,CC,GG,miua,para,sigma,c0,0,src);
%             end
%         end
%       %%%%   -------------------------------------------

        tcal = tcal + 1; msg = ['tcal = ',num2str(tcal)]; disp(msg);
        
%         %%%Error Calculations ----------------------------------------------          
%         if counter == 1
%             UU_fine(:,:,tcal) = UU; 
%             U_fine(:,:,tcal) = u; 
%             V_fine(:,:,tcal) = v; 
%             p_fine(:,:,tcal) = vn;
%         elseif counter == 2
%             maxs1 = max(max(abs(UU_fine(1:2:nfine, 1:2:mfine,tcal)-UU)));
%             maxp1 = max(max(abs(p_fine(1:2:nfine, 1:2:mfine,tcal)-vn)));
%             maxu1 = max(max(sqrt((U_fine(1:2:nfine, 1:2:mfine,tcal)-u).^2+(V_fine(1:2:nfine, 1:2:mfine,tcal)-v).^2)));
%             Smax_1 = max(Smax_1,maxs1); umax_1 = max(umax_1,maxu1); pmax_1 = max(pmax_1,maxp1);
%             l2_s1 = norm(UU_fine(1:2:nfine, 1:2:mfine,tcal)-UU,'fro');
%             l2_p1 = norm(p_fine(1:2:nfine, 1:2:mfine,tcal)-vn,'fro');
%             l2_u1 = sqrt(norm(U_fine(1:2:nfine, 1:2:mfine,tcal)-u,'fro').^2+ norm(U_fine(1:2:nfine, 1:2:mfine,tcal)-u,'fro').^2);
%             S2_c1= S2_c1+l2_s1^2; p2_c1= p2_c1+l2_p1^2; u2_c1= u2_c1+l2_u1^2;
%         elseif counter == 3
%             maxs2 = max(max(abs(UU_fine(1:4:nfine, 1:4:mfine,tcal)-UU)));
%             maxp2 = max(max(abs(p_fine(1:4:nfine, 1:4:mfine,tcal)-vn)));
%             maxu2 = max(max(sqrt((U_fine(1:4:nfine, 1:4:mfine,tcal)-u).^2+(V_fine(1:4:nfine, 1:4:mfine,tcal)-v).^2)));
%             Smax_2 = max(Smax_2,maxs2); umax_2 = max(umax_2,maxu2); pmax_2 = max(pmax_2,maxp2);
%             l2_s2 = norm(UU_fine(1:4:nfine, 1:4:mfine,tcal)-UU,'fro');
%             l2_p2 = norm(p_fine(1:4:nfine, 1:4:mfine,tcal)-vn,'fro');
%             l2_u2 = sqrt(norm(U_fine(1:4:nfine, 1:4:mfine,tcal)-u,'fro').^2+ norm(U_fine(1:4:nfine, 1:4:mfine,tcal)-u,'fro').^2);
%             S2_c2= S2_c2+l2_s2^2; p2_c2= p2_c2+l2_p2^2; u2_c2= u2_c2+l2_u2^2;
%         else
%             maxs3 = max(max(abs(UU_fine(1:8:nfine, 1:8:mfine,tcal)-UU)));
%             maxp3 = max(max(abs(p_fine(1:8:nfine, 1:8:mfine,tcal)-vn)));
%             maxu3 = max(max(sqrt((U_fine(1:8:nfine, 1:8:mfine,tcal)-u).^2+(V_fine(1:8:nfine, 1:8:mfine,tcal)-v).^2)));
%             Smax_3 = max(Smax_3,maxs3); umax_3 = max(umax_3,maxu3); pmax_3 = max(pmax_3,maxp3);
%             l2_s3 = norm(UU_fine(1:8:nfine, 1:8:mfine,tcal)-UU,'fro');
%             l2_p3 = norm(p_fine(1:8:nfine, 1:8:mfine,tcal)-vn,'fro');
%             l2_u3 = sqrt(norm(U_fine(1:8:nfine, 1:8:mfine,tcal)-u,'fro').^2+ norm(U_fine(1:8:nfine, 1:8:mfine,tcal)-u,'fro').^2);
%             S2_c3= S2_c3+l2_s3^2; p2_c3= p2_c3+l2_p3^2; u2_c3= u2_c3+l2_u3^2;
%         end
%         %%%%%%---------------------------------------------------------
        
        %%%% Cumulative Oil recovered calculations------------------------
        if tcal == 1
            COC(counter,tcal) = OC;
        else
            COC(counter,tcal) = COC(counter,tcal-1)+OC;
        end      
        CROIP(counter,tcal) = ROIP;
        %%%%% -----------------------------------------------------------
        
%         u_old = u;
%         v_old = v;

% %%%%%%%%%% Get Visualization of saturation/concentration profiles ---------
            figure(1)
            contourf(x,y,UU,5); shading flat; view(0,90); % view for contourf plots
            %contourf(x,y,UU,[0.24 0.46 0.7 0.95]); shading flat; view(0,90); % view for contourf plots
            %surf(x,y,UU); shading interp; view(105,50); % view for surf plots
            caxis manual; caxis([0.21 (1-sor0)]);colorbar; colormap jet;
            drawnow  %for live view of figures during simulation
%             fig1 = figure(1); fig1.OuterPosition = [0 0 900 1200];
%             
            %%%% Make a water saturation movie 
%             if counter == 1
%                 mov1 = getframe(gcf); writeVideo(vid1,mov1);
%             elseif counter == 2
%                 mov2 = getframe(gcf); writeVideo(vid2,mov2);
%             else
%                 mov3 = getframe(gcf); writeVideo(vid3,mov3);
%             end

%             if rem(tcal,5)==0
%                 export_fig(sprintf('./Sim6/Run%d/ST%dK%dx%d.pdf', counter, tcal, size(KK,1),size(KK,2)),'-opengl');
%             end
%
%
            if c0 ~=0
                figure(2)
                h2=contourf(x,y,CC);
                caxis manual; caxis([0 c0+0.02*c0]); view(0,90); colorbar; colormap spring; shading flat;
                drawnow  %for live view of figures during simulation
% %                 fig2 = figure(2); fig2.OuterPosition = [0 0 1080 1400]; 
% %                 
% % %                 if rem(tcal,50)==0
% % %                     export_fig(sprintf('./Sim5/Run%d/CT%dK%dx%d.pdf', counter, tcal, size(KK,1),size(KK,2)),'-opengl');
% % %                 end
% % 
            end

            if g0~=0
                figure(3)
                h3=contourf(x,y,GG,5);
                caxis manual; caxis([0 g0+0.02*g0]); view(0,90); colorbar; colormap cool; shading flat;
                drawnow  %for live view of figures during simulation
%                 %fig3 = figure(3); fig3.OuterPosition = [0 0 1080 1000];
%                 
% %                 if rem(tcal,100)==0
% %                     export_fig(sprintf('./Sim5/Run%d/GT%dK%dx%d.pdf',counter, tcal, size(KK,1),size(KK,2)),'-opengl');
% %                 end
            end
           
%%%%%%%%%%% ---------------------------------------------------------------
    end
    %%%% Time simulation loop ends here with water breakthrough or similar
    %%%% condition --------------------------------------------------------
    
    
%%%%%% Plots of oil recovery, area swept etc. ---------------------------    
    figure(20)
    scatter(1:tcal,COC(counter,1:tcal)); title('Cumulative Oil Recovered');
    xlabel('Time');ylabel('Volume of Oil Recovered');
%     export_fig(sprintf('./Sim2/COR%dG%0.3fK%dx%d.pdf', counter, g0 , size(KK,1),size(KK,2)));
    figure(21)
    scatter(1:tcal,CROIP(counter,1:tcal)); title('Residual Oil in place(ROIP)');
    xlabel('Time');ylabel('Volume of residual oil');
%     export_fig(sprintf('./Sim2/ROIP%dG%0.3fK%dx%d.pdf', counter, g0 ,size(KK,1),size(KK,2)));
    
     time(counter) =toc; iterations(counter) = tcal;
%     if counter == 1 
% 	time1 = time(1); iterations1 = tcal; COC1 = COC(1,1:tcal); CROIP1= CROIP(1,1:tcal);
% 	save ./Sim2/runtime-lenovosim2.mat time1 iterations1 COC1 CROIP1
% 	clear time1 iterations1 COC1 CROIP1
%     elseif counter == 2
% 	time2 = time(2); iterations2 = tcal; COC2 = COC(2,1:tcal); CROIP2= CROIP(2,1:tcal);
% 	save ./Sim2/runtime-lenovosim2.mat time2 iterations2 COC2 CROIP2 -append 
% 	clear time2 iterations2 COC2 CROIP2
%     else
% 	time3 = time(3); iterations3 = tcal; COC3 = COC(3,1:tcal); CROIP3= CROIP(3,1:tcal);
% 	save ./Sim2/runtime-lenovosim2.mat time3 iterations3 COC3 CROIP3 -append
% 	clear time3 iterations3 COC3 CROIP3
%     end 
%%%% ---------------------------------------------------------------------

% %%% Error calculations  --------------------------------------------------
%     if counter == 1
%         UU_fine(:,:,tcal:nFrames) = []; U_fine(:,:,tcal:nFrames) = []; p_fine(:,:,tcal:nFrames) = [];
%     elseif counter == 2
%         S2_c1 = sqrt(S2_c1*dt*para.box.dx*para.box.dy); l2_s1 = l2_s1*sqrt(para.box.dx*para.box.dy);
%         p2_c1 = sqrt(p2_c1*dt*para.box.dx*para.box.dy); l2_p1 = l2_p1*sqrt(para.box.dx*para.box.dy);
%         u2_c1 = sqrt(u2_c1*dt*para.box.dx*para.box.dy); l2_u1 = l2_u1*sqrt(para.box.dx*para.box.dy);
%     elseif counter == 3
%         S2_c2 = sqrt(S2_c2*dt*para.box.dx*para.box.dy); l2_s2 = l2_s2*sqrt(para.box.dx*para.box.dy);
%         p2_c2 = sqrt(p2_c2*dt*para.box.dx*para.box.dy); l2_p2 = l2_p2*sqrt(para.box.dx*para.box.dy);
%         u2_c2 = sqrt(u2_c2*dt*para.box.dx*para.box.dy); l2_u2 = l2_u2*sqrt(para.box.dx*para.box.dy);
%     else 
%         S2_c3 = sqrt(S2_c3*dt*para.box.dx*para.box.dy); l2_s3 = l2_s3*sqrt(para.box.dx*para.box.dy);
%         p2_c3 = sqrt(p2_c3*dt*para.box.dx*para.box.dy); l2_p3 = l2_p3*sqrt(para.box.dx*para.box.dy);
%         u2_c3 = sqrt(u2_c3*dt*para.box.dx*para.box.dy); l2_u3 = l2_u3*sqrt(para.box.dx*para.box.dy);
%     end
% %%%% ---------------------------------------------------------------------------------
    
% %%% Closing the video structure %%%%%%%%%%%%%%
% if counter == 1
%     close(vid1);
% elseif counter == 2
%     close(vid2);
% elseif counter == 3
%     close(vid3);
% end

end
% %%% If variables need to be stored ---------------------------------------
% save ./runtime-adasim.mat time iterations Smax_1 Smax_2 Smax_3 umax_1 pmax_1 umax_2 pmax_2 umax_3 pmax_3 S2_c1 S2_c2 S2_c3 
% save ./runtime-adasim.mat p2_c1 p2_c2 p2_c3 u2_c1 u2_c2 u2_c3 l2_s1 l2_s2 l2_s3 l2_p1 l2_p2 l2_p3 l2_u1 l2_u2 l2_u3 -append
% %%%% ------------------------------

% %%%% Profiling for code improvement
% % profile viewer
% % p=profile('info');
% % profsave(p,'profile_results')

disp('Simulation complete')
% exit   %Only use for Ada implementation

        