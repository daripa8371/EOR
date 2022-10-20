function [swr,sor] = compres(sigma,u,v,miua)
%function to compute residual saturations as a function of surfactant
%concentration via capillary number variation. Hence it varies with change
%in surf concentration and velocity evolution. Must be recomputed at every
%time step.

global miuo swr0 sor0;

%m = size(u,2); n = size(u,1);
% define critical capillary numbers 
% ie $$N_c $$ at which $$s_{ro}$$ and $$ s_{ra}$$ begin to decrease
Nco0 = 1.44*10^(-4);  %% Values from Amafuele Handy 1982
Nca0 = 1.44*10^(-4);  %%% these two do not have to be the same

% compute capillary number
nca = sqrt(u^2+v^2).*miua./sigma; nco = sqrt(u^2+v^2)*miuo./sigma;
Nca = norm(nca); % compute 2-norm of Nc matrix ie largest singular value
Nco = norm(nco);
% define residual saturations as functions of capillary numbers
if Nco< Nco0
    sor = sor0;
elseif Nco >= Nco0
    sor = sor0*(Nco0/Nco)^(0.5213);
end

if Nca < Nca0
    swr = swr0;
elseif Nca>= Nca0
    swr = swr0*(Nca0/Nca)^(0.1534);
end

end