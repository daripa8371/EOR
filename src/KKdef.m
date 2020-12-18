function [] = KKdef(~,para,x,y,flag)

global KK

switch flag
    case 1        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Defining a homogeneous permeability field
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        KK = 100*ones(para.box.n+1,para.box.m+1);
        
        %-------------------------------------------
        
    case 2
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Defining a continuous heterogeneous permeability function
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        KK = 1000*sin(4*x) + 1000*cos(5*y) + 100*exp((x-.5)^2+(y-.5)^2);
        
        %KK = 
        
        
        %     %-------------------------------------------------------------
        
    case 3
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Defining a permeability function with an impermeable block
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        bnx = para.box.m+1; bny = para.box.n+1;
        KK = 3000*ones(bny,bnx);
        KK(floor(bny/2)-floor(bny/8):floor(bny/2)+floor(bny/8),floor(bnx/2)-floor(bnx/8):floor(bnx/2)+floor(bnx/8)) = 3;
        
    case 4
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Defining a permeability function with two impermeable blocks
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        bnx = para.box.m+1; bny = para.box.n+1;
        KK = 3000*ones(bny,bnx);
        KK(floor(3*bny/4)-floor(bny/12):floor(3*bny/4)+floor(bny/12),floor(2*bnx/3)-floor(bnx/12):floor(2*bnx/3)+floor(bnx/12)) = 3;  
        KK(floor(bny/3)-floor(bny/10):floor(bny/3)+floor(bny/10),floor(bnx/3)-floor(bnx/10):floor(bnx/3)+floor(bnx/10)) = 3;    
        
    case 5     
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Simulating with Upper Ness permeability from SPE10
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        load('KK30Ness.mat'); disp('Upper Ness formation permeability loaded')
    case 6
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Simulating with Tabert permeability from SPE10
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        load('KK30Tabert.mat'); disp('Tarbert formation permeability loaded')
        %---------------------------------------------------------------
    case 7
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Simulating with new Tarbert permeability from SPE10
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        load('KK30Tarbertnew.mat'); disp('Tarbert formation permeability loaded')
        %---------------------------------------------------------------

    case 8
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Simulating with new Upper Ness permeability from SPE10
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        load('KK30Nessnew2.mat'); disp('Upper Ness formation permeability loaded')
        %---------------------------------------------------------------
    case 9
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Permeability function from  Ferreira, Pena, Romanazzi 2015
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %bn = sizeofgrid(counter)+1;
        KK = 100*0.5*( 0.5*(1-10^(-7))*(sin(6*pi*cos(x)).*cos(4*pi*sin(3*y))-1)+1);
        
        %---------------------------------------------------------------
    case 10
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Defining a permeability function with impermeable block channels
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        bnx = para.box.m+1; bny = para.box.n+1;
        KK = 3000*ones(bny,bnx);
        KK(floor(bny/9):floor(bny/2),floor(bnx/6):floor(bnx/6)+floor(bnx/20)) = 3;  
        KK(floor(bny/3):bny,floor(bnx/2):floor(bnx/2)+floor(bnx/15)) = 3; 
        KK(floor(2*bny/3):floor(5*bny/6),floor(5*bnx/6):floor(5*bnx/6)+floor(bnx/20)) = 3;
    case 11
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Permeability function from  Chueh, Djilali, Bangerth 2013
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        KK = max(100*exp(-((y-0.5-0.1*sin(10*(x)))/0.1).^2),0.01);    
        
    case 12
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Simulating with new multiscale Gaussian permeability 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        load('K32new.mat'); disp('Multiscale permeability field loaded')
        %---------------------------------------------------------------
        
    case 13
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Grid Orientation permeability 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        KK = 3.72e-3 * (1 + 0.01*cos(50*pi*(x-0.5)).*cos(50*pi*(y-0.5)));
        
    case 14
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Grid Orientation permeability
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        KK = 3.72e-3*ones(size(y,1),size(x,2));
        
      %%%%%%%%%%%% Plotting the permeability field %%%%%%%%%%%%%%%%%%
%     figure(100)
%     surf(x,y,KK);
% %     export_fig(sprintf('KK%dx%d.pdf',para.box.m+1, para.box.n+1),'-opengl');
% %    
%     figure(101)
%     pcolor(KK); shading interp; colormap jet; colorbar; title('Ferreira, Pena, Romanazzi (2015)'); drawnow
%     export_fig(sprintf('KKFPR2015-raster%dx%d.pdf',para.box.m+1, para.box.n+1),'-opengl');
%      pause
%     %---------------------------------------------------------------


end    



%     %---------------------------------------------------------------
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % Defining a random matrix for heterogeneous coefficient
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %
%     if counter == 1
%         load('KK16.mat');
%     elseif counter == 2
%         load('KK16.mat');
%         KK1 = interp2(KK,1); KK = KK1; clear KK1;
%     else 
%         load('KK16.mat')
%         KK2 = interp2(KK,2); KK = KK2; clear KK2;
%     end
%     %---------------------------------------------------------------

    

    