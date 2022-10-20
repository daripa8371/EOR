function lambda=compmob(s,miua,c,sor,swr,flag,surf)
%%% function to compute mobility 
% sor = residual oil saturation at IFT sigma (matrix)
% swr = residual water saturation at IFT sigma (matrix)
% sor0 = residual oil saturation at IFT sigma0 (constant)
% swr0 = residual water saturation at IFT sigma0 (constant)
% flag = denotes which phase 0=oleic and 1=aqueous
% surf = flag that denotes presence (1) or absence (0) of surfactant

global sor0 swr0 miuo;
% 
% miua = compvis(c);
%%
if surf == 0
%normalized saturations of water and oil at IFT sigma0
nsw0 = (s-swr0)/(1-swr0);
nso0 = (s-swr0)/(1-swr0-sor0); 

%Corey type relative permeability in absence of surfactant
krw0 = nsw0.^3.5;
kro0 = ((1-nso0).^2).*(1-nso0.^(1.5));

if flag == 0 %% oil
    lambda = kro0./miuo;
elseif flag == 1 %% water
    lambda = krw0./miua;
end
%%
elseif surf == 1
%normalized saturations of water and oil at IFT sigma
nsw = (s-swr)/(1-swr);
nso = (s-swr)/(1-swr-sor);

%rel perm in presence of surfactant
krw = nsw.*(2.5*swr.*(nsw.^2-1)+1);
kro = (1-nso).*(1-5*sor.*nso);


if flag == 0 %% oil
    lambda = kro./miuo;
elseif flag == 1 %% water
    lambda = krw./miua;
end
end
%%