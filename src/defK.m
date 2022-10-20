clc 
clear all
% matlab structure array for defining mesh and grid points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
para.box.left=0;
para.box.right=1;
para.box.bottom=0;
para.box.top=1;
para.box.m=23;     % x grid points is 0 to para.box.m
para.box.n=23;     % y grid points is 0 to para.box.n
para.box.dx = (para.box.right-para.box.left)/para.box.m;
para.box.dy = (para.box.top-para.box.bottom)/para.box.n;
[x,y] = meshgrid(para.box.left:para.box.dx:para.box.right,para.box.bottom:para.box.dy:para.box.top);

K=randn(para.box.n+1, para.box.m+1);
K=exp(K);
% KK =ones(para.box.n+1, para.box.m+1);  %%% For a constant permeability
figure(10)
surf(x,y,K);
view(-30,30);
figure(11)
contourf(x,y,K);
% [C,h] = contour(x,y,K,10);
% get(h,'LevelStep')
% set(h,'ShowText','on','TextStep',get(h,'LevelStep')*2)
%        contour(x,y,K);
%        colorbar
