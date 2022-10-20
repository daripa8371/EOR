function [] = KKdef(counter,sizeofgrid,x,y,flag)

global KK

switch flag
    case 1        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Defining a homogeneous permeability field
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        Kmax = 100; %JCP 2017 Daripa and Dutta used 1000
        KK = Kmax*ones(sizeofgrid(counter)+1);
        
        %-------------------------------------------
        
    case 2
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Defining a continuous heterogeneous permeability function
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %        KK = 1000*sin(4*x) + 1000*cos(5*y) + 100*exp((x-.5)^2+(y-.5)^2);
        
          Kmax = 100; %JCP 2017 Daripa and Dutta used 50
          KK = Kmax*( 0.5*(1-10^(-7))*(sin(6*pi*cos(x)).*cos(4*pi*sin(3*y))-1)+1);
        %     %-------------------------------------------------------------
        
    case 3
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Defining a permeability function with an impermeable block
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        bn = sizeofgrid(counter)+1;
        KK = 3000*ones(bn);
        KK(floor(bn/2)-floor(bn/8):floor(bn/2)+floor(bn/8),floor(bn/2)-floor(bn/8):floor(bn/2)+floor(bn/8)) = 3;
        
    case 4
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Defining a permeability function with an impermeable block
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        bn = sizeofgrid(counter)+1;
        KK = 3000*ones(bn);
        KK(floor(3*bn/4)-floor(bn/12):floor(3*bn/4)+floor(bn/12),floor(2*bn/3)-floor(bn/12):floor(2*bn/3)+floor(bn/12)) = 3;  
        KK(floor(bn/3)-floor(bn/10):floor(bn/3)+floor(bn/10),floor(bn/3)-floor(bn/10):floor(bn/3)+floor(bn/10)) = 3;    
        
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

end    
%     %%%%%%%%%%%% Plotting the permeability field %%%%%%%%%%%%%%%%%%
%     figure(100)
%     surf(x,y,KK);
%     export_fig(sprintf('KK%dx%d.pdf',para.box.m+1, para.box.n+1),'-opengl');
%     
%     pause
%     %---------------------------------------------------------------



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

    

    