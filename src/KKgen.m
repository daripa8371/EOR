clear variables
clc

para.box.m=15;     % x grid points is 0 to para.box.m
para.box.n=15;     % y grid points is 0 to para.box.n
para.box.left=0;
para.box.right=1;
para.box.bottom=0;
para.box.top=1;
para.box.dx = (para.box.right-para.box.left)/para.box.m;
para.box.dy = (para.box.top-para.box.bottom)/para.box.n;
[x,y] = meshgrid(para.box.left:para.box.dx:para.box.right,para.box.bottom:para.box.dy:para.box.top);

 KK=randn(para.box.n+1, para.box.m+1);
    KK=exp(KK);
    figure(10)
    surf(x,y,KK);