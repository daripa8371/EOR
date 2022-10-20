%% Code to implement shear effects of polymer on the EOR simulation
%% Rohit Mishra
%% 04/21/2019
%% Todd-Longstaff mixing model

function [mu_w_eff,m_mu] = TLmixing(miua,c)

global miuw miup_array c0 miuo beta1
omega = 0;     

for ii=1:size(miua,1)
    for jj=1:size(miua,2)
        if (c(ii,jj)>0) 
            
            mu_p_eff(ii,jj) = (miua(ii,jj)^omega) * (miup_array(ii,jj)^(1-omega));
            mu_w_e(ii,jj) = (miua(ii,jj)^omega) * (miuw^(1-omega));
            C_bar(ii,jj) = c(ii,jj)/c0;
            mu_w_eff (ii,jj) = mu_w_e (ii,jj)* mu_p_eff (ii,jj)/ ...
                ( (1-C_bar(ii,jj)) * mu_p_eff(ii,jj)...
                + C_bar(ii,jj) * mu_w_e(ii,jj)  );
            
            %% Calculating shear multipliers
            m_mu(ii,jj) = (1+beta1*c(ii,jj));
        else
            mu_w_eff(ii,jj) = miua(ii,jj);
            m_mu(ii,jj) = 1;
        end
        
    end
end
end

