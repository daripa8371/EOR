function [c0]=set_c0_rev_fing(para,phi,c_0)
% function to initialize a viscous profile in the displacing phase 
% to test reverse fingering phenomenon

global init_front_hs

m = para.box.m;
n = para.box.n;

dx = para.box.dx;
dy = para.box.dy;

left = para.box.left;
bottom = para.box.bottom;

E = (phi<0);
[ii,jj] = meshgrid(1:m+1,1:n+1);

x = left+(ii-1)*dx; 
y = bottom+(jj-1)*dy;

C0 = (0.3-c_0).*x/init_front_hs + c_0;
c0 = E.*C0;
end