function miua=compvis(c)
% function to compute viscosity of injected 
% displacing phase containing polymer
% miuo = Displaced phase viscosity, c = polymer concentration

global miuw miuo c0

beta1 = 0.5;

n=size(c,1); m=size(c,2);
if c0 == 0
    miua = miuw*ones(n,m);
else
    miua = miuw*(1+beta1*c); %miuo*(0.5+c);%
end