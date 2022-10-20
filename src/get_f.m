function [out,out_]=get_f(para,beta,beta_,phi)
dx=para.box.dx;
dy=para.box.dy;
m=para.box.m;
n=para.box.n;
f=para.f;
f_=para.f_;

[px,py]=get_gra(beta,para,phi);
% [px_,py_]=gradient(beta_,dx,dy);
% [px_,py_]=get_gra_con(beta_,para);
% [px,py,px_,py_]=get_beta_xy(para,phi);
out=zeros(n+1,m+1);
for i=1:m+1
    for j=1:n+1
        x=para.box.left+(i-1)*dx;
        y=para.box.bottom+(j-1)*dy;
        if phi_func(x,y,para,phi)>10^(-10)||abs(phi_func(x,y,para,phi))<10^(-10)
            out(j,i)=f(x,y,beta_func(x,y,para,beta),[px(j,i);py(j,i)]);
        else
            out(j,i)=f_(x,y,beta_func(x,y,para,beta),[px(j,i);py(j,i)]);
        end
        out_(j,i)=f_(x,y,beta_func(x,y,para,beta_),[px(j,i);py(j,i)]);
    end
end