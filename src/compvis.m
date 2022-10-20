
function [miua,gamma_dot]=compvis(c,U,V,X,Y)
% function to compute viscosity of injected
% displacing phase containing polymer
% miuo = Displaced phase viscosity, c = polymer concentration
%

global viscosityFlag
global miuw miuo beta1 c0 miup_array c0_array polymerType

gamma_dot =zeros;
switch viscosityFlag
    case 1
        
        
        
        %% Base code
        
        
        
        n=size(c,1); m=size(c,2);
        if c0 == 0
            miua = miuw*ones(n,m);
        else
            miua = miuw*(1+beta1*c); %miuo*(0.5+c);%
        end
        
    case 2
        
        % Sourav's implementation
        n=size(c,1);
        if c0 == 0
            miua = miuw*ones(n);
        else
            miua = miuo*(0.5+c);
        end
        
        
        
        
    case 3        
        % Modified Nejat-Basagaoglu
        rho_water = 1000;
        rho_xanthane = 1500;
      
        w1=rho_xanthane * c;
        w2=rho_water*(1-c);
        wppm = ((w1)./ (w1 +w2))*1e6;
        w10 = rho_xanthane * c0_array;
        w20=rho_water*(1-c0_array);
        wppm0 = ((w10)./ (w10 +w20))*1e6;
        
%         [ 4.86265534 -0.41570227] (Schizo – n)
%         [0.03647214 1.32175949] (Schizo – epsilon)
%                 
%         [ 3.05428284 -0.27294817] (Xanthane – n)
%         [1.15410398e-04 2.04937780e+00] (Xanthane – epsilon)
    if(polymerType)
        
       epsCoeff = [0.03647214 1.32175949];
       nCoeff =  [ 4.86265534 -0.41570227];
       
    else
        

       epsCoeff = [1.15410398e-04; 2.04937780e+00];
       nCoeff =  [ 3.05428284; -0.27294817];

    end
     




       epsilon0 = epsCoeff(1)*wppm0.^epsCoeff(2);
        power_n0 = min(nCoeff(1)*wppm0.^nCoeff(2),1);
       
        epsilon = epsCoeff(1)*wppm.^epsCoeff(2);
        power_n = min(nCoeff(1)*wppm.^nCoeff(2),1);
        
%         epsilon0(:) = 250.*ones;
%         power_n0(:) = 0.5.*ones;
%         epsilon(:) = 250.*ones;
%         power_n(:) = 0.5.*ones;
        

         n=size(c,1); m=size(c,2);
        miua = miuw*ones(n,m);
     
        a1=(divergence(X,V));
        a2=(divergence(Y,U));
        a3=(divergence(X,U));
        a4=(divergence(Y,V));
        pi_D = abs(-0.25*((a1+a2).^2)+a3*a4);
        
        for ii=1:size(miua,1)
            for jj=1:size(miua,2)
                if (c(ii,jj)>0)
                    
                    gamma_dot(ii,jj) = 2*sqrt(pi_D(ii,jj));
                    if gamma_dot(ii,jj) == 0
                        
                       
                    else
                        miup_array(ii,jj) = epsilon0(ii,jj) * (gamma_dot(ii,jj)^(power_n0(ii,jj)-1));
                        
                        miua(ii,jj) = epsilon(ii,jj) * (gamma_dot(ii,jj)^(power_n(ii,jj)-1));
                        if (miua(ii,jj) < miuw)
                            miua(ii,jj) = miuw;
                        end
                         if (miua(ii,jj) > 100)
                            miua(ii,jj) = 100;
                         end
                         if (miup_array(ii,jj) < miuw)
                            miup_array(ii,jj) = miuw;
                        end
                        if (miup_array(ii,jj) > 100)
                            miup_array(ii,jj) = 100;
                        end
                    end
                end
            end
        end
        
     gammaMax = max(max(gamma_dot));   
end
end



