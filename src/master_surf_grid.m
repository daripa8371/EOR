%% Main program for a full flooding simulation
% Grid Sizes used 15x15, 20x20, 30x30, 40x40 
clc
clear variables;
clear global;
format long
profile clear
profile on



% declaring global variables to reduce computational time in passing
% them through function calls
global c0_array miuo miuw miup swr0 sor0 dt KK ...
    s0 c0 g0 beta1 viscosityFlag shearFlag miup_array ...
         polymerType  ;

% [obs,obs_header,obs_title]=read_eas('syn_perm.out');

%%%% matlab structure array for defining mesh and grid points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%
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
nsim =1;
sog=29;
sizeofgrid = [sog sog sog sog];
c0iter = [0.001 0 0 0.001 ]; %0.0002 (300 wppm) 0.0006 (900 wpppm) 0.001 (1500 wpppm)
g0iter = [0 0 0 0]; %[0.01 0.01 0.01];
k = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% preallocating variables that will store cpu time taken for 
%%% each run and number of iterations/time-steps of each run
time = zeros(1,3); iterations = zeros(1,3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for counter = 1:nsim %  3:4    
    tic
    
    para.box.m  = sizeofgrid(counter);     % x grid points is 0 to para.box.m
    para.box.n  = sizeofgrid(counter);     % y grid points is 0 to para.box.n
    para.box.dx = (para.box.right-para.box.left)/para.box.m;
    para.box.dy = (para.box.top-para.box.bottom)/para.box.n;
    [x,y] = meshgrid(para.box.left:para.box.dx:para.box.right,para.box.bottom:para.box.dy:para.box.top);
 
    %%% matrix containing the evaluation of the distance based level set
    % function at each of the grid points
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    phi_test=get_phi_test(para); 
        
    %Defining right hand side of elliptic system - source terms
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    f=zeros(para.box.n+1,para.box.m+1);
    MFW=zeros;
    iterX_save = zeros;
    
     src = 120000; %50000; %5000; %2; %200 for QFS    % magnitude of mass flow rate at source (1.2 for heterogeneous)1.2 to 30 for shear cases
%     rate = 0; 
    %% Rectilinear propagation
    f(:,1)=src;                           % intensity of injection well = src
    f(:,para.box.m+1)=-src;    % intensity of production well = -src

%     %% Quarter-five spot
%   f(1,1)=src;                           % intensity of injection well = src
%   f(para.box.n+1,para.box.m+1)=-src;    % intensity of production well = -src
    %-------------------------------------------------------------------------
    
    
    %%%%%%%   Defining permeability matrix    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % flag values:  1 = homogeneous with magnitude 1000
    %               2 = continuous heterogeneous function
    %               3 = impermeable block inclusion at the center
    %               4 = impermeable block inclusion off the center
    %               5 = Upper Ness from SPE10 model sections
    %               6 = Tabert from SPE10 model sections 
    
    permeabilityFlag = 1;
    KKdef(counter,sizeofgrid,x,y,permeabilityFlag);

%     KK = 100*KK;
    %-------------------------------------------------------------------
    
    
    
    %%%%%%%%%%% Preparing automatic runs  %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % initialize RP variables----------------------------------------------
    % UU=zeros(para.box.n+1,para.box.m+1);
    s0=0.79;  % initial residual water saturation inside the reservoir = 1-s0
    c0=c0iter(counter);%0.1;
    g0=g0iter(counter);%0.005;
    c0_array = c0*ones(sog+1,sog+1);
    [UU,CC,GG]=s0c0(para,phi_test,s0,c0,g0,1);
    interface = zeros(30,1);
    %%%%%
    miuo=10; %2e-3;%12.6;  %%% 0.95
    miuw=1.26; %8.9e-4;%1.26;%1.26; %0.0430;   %%% 0.095   % Effective aq. viscosity = miuo*(0.5+c)
    beta1= 15000; %2;
    miup = miuw*(1+beta1*c0); %% Assumption
    miup_array = miup*ones(sog+1,sog+1);
    %alpha=0.6;  %%% 0.6
    % define initial residual saturations before critical capillary number 
    swr0 = 0.1; % actual initial residual saturation is higher than this value in practice
    sor0 = 0.3;
    
    viscosityFlag = 3; %Flag to enable different viscosity models
    % 3 = Modified Nejat-Basagaoglu model (set src to 5000)
    % 2 = Sourav's model
    % 1 = Original model from JCP paper
    
    polymerType = 1;
    % 0 = Xanthane
    % 1 = Schizophyllan
    
    TLmixingFlag = 0; %Flag to enable Todd-Longstaff mixing parameter
    % 1 = TL mixing model ON
    % 0 = TL mixing model OFF
   
    shearFlag = 0; %Flag to enable shear effects with ad-hoc model from Eclipse code
    % 1 = Shear effects ON
    % 0 = Shear effects OFF
    

    t=0; tcal = 0;
    %cfl=0.2;
    u = zeros(para.box.n+1,para.box.m+1); v = u;
    
    %     u_old = zeros(para.box.n+1,para.box.m+1); v_old=u_old; u=u_old; v=v_old;
    CFL = 1;
    dt = CFL * para.box.dx / src;     %* 100; only for quarter five spot
    % dt = 1/25;
    tstop=500;
    COC = zeros(nsim,2000); %floor(tstop/dt));
    ProdRate = zeros(nsim,floor(tstop/dt));
    CROIP = zeros(nsim,floor(tstop/dt));
    miuaSave = zeros;
    shearSave = zeros;
    concSave = zeros;
    src_total = 0;
    sumUU = 0;
    while (t<tstop && UU(para.box.n+1,para.box.m+1) <= 0.70)

      %      src = src + rate*tcal*0.22941; %200 for QFS    % magnitude of mass flow rate at source (1.2 for heterogeneous)1.2 to 30 for shear cases
         src_total = src_total + src;
     
        t=t+dt;
        innerIter = 0;
        epsilon = 10;
        %Solve elliptic equation.................................
%        while(epsilon > 1e-4) 
        % IFT as a function of surf concentration
        % $$ \sigma = \frac{10.001}{\Gamma +1} - 0.001 $$
        sigma = 10.001./(GG+1)-0.001;
        
        % aq soln viscosity as a function of polymer $$ \mu_a =
        % \mu_o(0.5+c) $$
        [miua,shear] = compvis(CC,u,v,x,y);
       if (viscosityFlag == 1)
           miup = max(miua(1,:));
        miup_array = miup*ones(sog+1,sog+1);
       end
        if (TLmixingFlag)
        [mu_w_eff,m_mu] = TLmixing(miua,CC);
        miua = mu_w_eff;
        end
        tSave = tcal+1;
        
[mShear, nShear] = size(shear);

        % recompute residual saturations $$ s_{ra}, s_{ro} $$
        [swr,sor] = compres(sigma,u,v,miua);
 
        % recompute mobilities ( with surfactants )
        lambda_a = compmob(UU,miua,CC,sor,swr,1,1);
        lambda_o = compmob(UU,miua,CC,sor,swr,0,1);
        
        lambda = lambda_a + lambda_o;
        
        beta=KK.*lambda;  %%% for polymer flow
        %     beta=KK.*((UU.^3)./(miuw)+((1-UU).^3)./miuo);  %%% for non-polymer
        
        [U,L] = setTri(para);
        gridsize = setGrid(para,U,L,beta);
        rh = setRightHand(para,f,U,L);
        A = setA(para,gridsize);
        A = spconvert(A);
        B = setB(para,gridsize,rh);
        uOld = u;
        vOld = v;
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
    epsilon = abs(max(u-uOld));
    if (innerIter > 1000)
        break;
    end
    innerIter = innerIter + 1;
    innerIterSave(tcal+1) = innerIter;
%        end
        
        %miup_array = compvis(c0_array,u,v,x,y);
        if (shearFlag)
       [miua] = shear_effects(u,v,miua,x,y,c);
        end
        
        
        %Solve transport problem....................................
        
        
%         [UU,CC,GG,OC,WC,ROIP]=nmmoc_surf_mod(u,v,UU,CC,GG,miua,para,sigma,c0,g0,src);
        
        [UU,CC,GG,OC,WC,ROIP]=nmmoc_surf_mod_neumann(u,v,UU,CC,GG,miua,para,sigma,c0,g0,src);

        tcal = tcal + 1; msg = ['tcal = ',num2str(tcal)]; disp(msg);
        % Cumulative Oil recovered
        if tcal == 1
            COC(counter,tcal) = OC;
        else
            COC(counter,tcal) = COC(counter,tcal-1)+OC;
        end
        ProdRate(counter,tcal) = OC/dt;
        CROIP(counter,tcal) = ROIP;
 

%%%%%%%%%% Get Visualization of saturation/concentration profiles ---------
           figure(1)
            miua_plot = miua/max(max(miua));
           contourf(x,y,miua); 
           % caxis([1.5,10])
            colorbar; colormap jet;
           drawnow  %for live view of figures during simulation
            

            
            figure(2)
            
           contourf(x,y,UU); %shading flat; view(0,90);
           colorbar; colormap jet;
           set(gca,'FontSize',18);

           title(tcal)
             [iterX_save(tcal), MFW(tcal)]=postProcess_fingerWidthQFS(UU);
%            [interface, MFW(tcal),iterX_save(tcal)]=postProcess_fingerWidth(UU);
%            if tcal>1
%            if(iterX_save(tcal) == iterX_save(tcal-1))
%                MFW(tcal) = max(MFW(tcal-1),MFW(tcal));
%            end
%            end
           %figure(3)
           %plot(interface,1:29);
            

            %% Code to make a movie 
            
            if counter == 1
                mov1(tcal) = getframe(gcf);
            elseif counter == 2
                mov2(tcal) = getframe(gcf);
            else
                mov3(tcal) = getframe(gcf);
            end
            
                        % pause


% saveLocation = strcat('./LogData/interface',string(tcal));
% 
            if rem(tcal,1)==0
                export_fig(sprintf('./Run%d/ST%dK%dx%d.jpeg', counter, tcal, size(KK,1),size(KK,2)),'-opengl');
%            save(saveLocation,'interface','-ascii');
            end

            
            if (rem(tcal,200) == 0)
                k = k+1;
                miuaTcal(:,:,k) = miua;
            end
%             if (tcal == 500 || tcal == 1000 || tcal == 1500)
%                 %hello
%                 miuanew=miua';
%                 tcal = tcal
%             end
%             
%%%%%%%%%%% ---------------------------------------------------------------
    end
    %%%% Time simulation loop ends here with water breakthrough or similar
    %%%% condition --------------------------------------------------------
    
    
%%%%%%% Plots of oil recovery, area swept etc. ---------------------------    
    figure(20)
    scatter(1:tcal,COC(counter,1:tcal)); title('Cumulative Oil Recovered');xlabel('Time');ylabel('Volume of Oil Recovered');
%    export_fig(sprintf('./Sim2/COR%dG%0.3fK%dx%d.pdf', counter, g0 , size(KK,1),size(KK,2)));
%     figure(21)
%     scatter(1:tcal,CROIP(counter,1:tcal)); title('Residual Oil in place(ROIP)');xlabel('Time');ylabel('Volume of residual oil');
%     export_fig(sprintf('./Sim2/ROIP%dG%0.3fK%dx%d.pdf', counter, g0 ,size(KK,1),size(KK,2)));
%     
%     time(counter) =toc; iterations(counter) = tcal;
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
%      end 
%%%% ---------------------------------------------------------------------

% ind = tcal+1:nFrames;
% if counter == 1
%     mov1(ind)=[];
% elseif counter == 2
%     mov2(ind) = [];
% elseif counter == 3
%     mov3(ind) = [];
% end

%size(mov)
end
% plot(iterX_save);
% save ./Sim2/runtime-lenovosim2.mat time iterations COC CROIP -append
% %%%% Profiling for code improvement
% % profile viewer
% % p=profile('info');
% % profsave(p,'profile_results')
% 
% %movie(M1);
% 
% movie2avi(mov1,'Block_inclusion.avi');
% movie2avi(mov2,'30Ness.avi');
% movie2avi(mov3,'30Tabert.avi');

% %% Viscosity plot
% [m1,n1] = size(concSave);
% for i=1:n1
%     for j = 1:n1
%         miu3D(i,j) = miua(j);
%     end
% end

        
        