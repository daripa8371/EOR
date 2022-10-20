function out=get_beta(para,v)


m = para.box.m;
n = para.box.n;

dx = para.box.dx;
dy = para.box.dy;

left = para.box.left;
bottom = para.box.bottom;

get_beta=zeros(n+1,m+1);
for ii=1:m+1
    for jj=1:n+1
           if v==1
            get_beta(jj,ii)=func_of_beta(left+(ii-1)*dx,bottom+(jj-1)*dy);
           end
           if v==-1
            get_beta(jj,ii)=func_of_beta_(left+(ii-1)*dx,bottom+(jj-1)*dy);
           end   
    end
end
out=get_beta;