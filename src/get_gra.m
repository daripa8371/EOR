function [px,py]=get_gra(vn,para)
  
m=para.box.m;
n=para.box.n;
dx=para.box.dx;
dy=para.box.dy;
px=zeros(n+1,m+1);
py=px;

for i=1:m+1
    for j=1:n+1
       
            if i~=1
                px(j,i)=(vn(j,i)-vn(j,i-1))/dx;
            end
            if i~=m+1
                px(j,i)=(vn(j,i+1)-vn(j,i))/dx;
            end
            if i~=1&&i~=m+1
                px(j,i)=(vn(j,i+1)-vn(j,i-1))/2/dx;
            end
            if j~=1
                py(j,i)=(vn(j,i)-vn(j-1,i))/dy;
            end
            if j~=n+1
                py(j,i)=(vn(j+1,i)-vn(j,i))/dy;
            end
            if j~=1&&j~=n+1
                py(j,i)=(vn(j+1,i)-vn(j-1,i))/2/dy;
            end
   
    end
end
                
  