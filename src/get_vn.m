function vn=get_vn(u,para)
  
m=para.box.m;
n=para.box.n;

vn=zeros(n+1,m+1);
for ii=1:m+1
    for jj=1:n+1
            vn(jj,ii)=u((jj-1)*(m+1)+ii);
    end
end
 