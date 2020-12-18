clear variables
clc

para.box.m=23;     % x grid points is 0 to para.box.m
para.box.n=23;     % y grid points is 0 to para.box.n
para.box.left=0;
para.box.right=1;
para.box.bottom=0;
para.box.top=1;
para.box.dx = (para.box.right-para.box.left)/para.box.m;
para.box.dy = (para.box.top-para.box.bottom)/para.box.n;
[x,y] = meshgrid(para.box.left:para.box.dx:para.box.right,para.box.bottom:para.box.dy:para.box.top);

 KK=randn(para.box.n+1, para.box.m+1);
    KK=34.4253*exp(KK);
    figure(10)
    surf(x,y,KK); shading interp;colormap jet; colorbar; title('Random permeability field')
    export_fig 'KK24x24.pdf' '-opengl'
    figure(11)
    R = log(KK/34.4253);
    pcolor(R); shading interp; colormap jet; colorbar; title('Gaussian field')
    export_fig 'Kraster24x24.pdf' '-opengl'