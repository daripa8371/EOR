function out=get_phi(para)


m = para.box.m;
n = para.box.n;

dx = para.box.dx;
dy = para.box.dy;

left = para.box.left;
bottom = para.box.bottom;

get_phi=zeros(n+1,m+1);
for ii=1:m+1
    for jj=1:n+1

            get_phi(jj,ii)=z_func(left+(ii-1)*dx,bottom+(jj-1)*dy);

    end
end
out=get_phi;