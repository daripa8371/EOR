function lambda=compmob(s,c,sor,swr,flag)
%%% function to compute mobility 
% sor = residual oil saturation at IFT sigma (matrix)
% swr = residual water saturation at IFT sigma (matrix)
% sor0 = residual oil saturation at IFT sigma0 (constant)
% swr0 = residual water saturation at IFT sigma0 (constant)
% flag = denotes which phase 0=oleic and 1=aqueous

global sor0 swr0 miuo g0 theta; %s0

miua = compvis(c);

%  if g0 == 0
%     %normalized saturations of water and oil at IFT sigma0
%     nsw0 = (s-swr0)/(1-swr0-sor0);
% %     nso0 = (s-swr0)/(1-swr0-sor0);
% 
%     %Corey type relative permeability in absence of surfactant
%     krw0 = nsw0.^((2+3*theta)/theta);
%     kro0 = ((1-nsw0).^2).*(1-nsw0.^((2+theta)/theta));
% %     kro0 = ((1-nso0).^2).*(1-nso0.^((2+theta)/theta));
% 
% %     nsw0 = (s-swr0)/(1-swr0-sor0);
% %     nso0 = (1-s-sor0)/(1-swr0-sor0);
% %     Mow = (1-s0-swr0)/(1-swr0-sor0);
% %     Omega = Mow+ nso0;    
% 
% %     krw0 = -999*ones(size(c,1),size(c,2));kro0=krw0;
% %     for j = 1:size(c,1)
% %         for i = 1:size(c,2)
% %             if s(j,i) <= Mow
% %                 krw0(j,i) = nsw0(j,i)^((2+3*theta)/theta);
% %                 kro0(j,i) = ((1-nso0(j,i))^2)*(1-nso0(j,i)^((2+theta)/theta));
% %             elseif s(j,i) > Mow
% %                 krw0(j,i) = nsw0(j,i)^2*(1+Mow^((2+theta)/theta) - Omega(j,i)^((2+theta)/theta));
% %                 kro0(j,i) = (1-nsw0(j,i))^2*(Omega(j,i)^((2+theta)/theta)-Mow^((2+theta)/theta));
% %             end
% %         end
% %     end
% %     fprintf('krw0(1,1) = %12.10f\n',krw0(1,1))
% %     fprintf('kro0(1,1) = %12.10f\n',kro0(1,1))  
%     
%     if flag == 0 %% oil
%         lambda = kro0./miuo;
%     elseif flag == 1 %% water
%         lambda = krw0./miua;
%     end
%     
%  elseif g0 ~= 0
%     %normalized saturations of water and oil at IFT sigma
%     nsw = (s-swr)/(1-swr);
%     nso = (s-swr)/(1-swr-sor);
    nsw = (s-swr)/(1-swr-sor);

%     
%     %rel perm in presence of surfactant
%     krw = nsw.*(2.5*swr.*(nsw.^2-1)+1);
%     kro = (1-nso).*(1-5*sor.*nso);  
    krw = nsw.^((2+3*theta)/theta);
    kro = ((1-nsw).^2).*(1-nsw.^((2+theta)/theta));

%     
    if flag == 0 %% oil
         lambda = kro./miuo;
    elseif flag == 1 %% water
        lambda = krw./miua;
    end

%     fprintf('krw(1,1) = %12.10f\n',krw(1,1))
%     fprintf('kro(1,1) = %12.10f\n',kro(1,1))
% end