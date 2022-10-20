clc;
clear;
c = 0.001;
rho_water = 1000;
rho_xanthane = 1500;
w1=rho_xanthane * c;
w2=rho_water*(1-c);
wppm1 = (w1)/ (w1 +w2)*1e6;
load xanthane_prop;
f1=fit(xanthane_prop(:,1),xanthane_prop(:,2),'poly2');
f2=fit(xanthane_prop(:,1),xanthane_prop(:,3),'poly2');


plot(f1,xanthane_prop(:,1),xanthane_prop(:,2));
