function out=get_new(u,para)
[nn,mm]=size(u);
sum=0;
for ii=1:nn
    sum=sum+u(ii);
end
aver=sum/nn;
out=u-aver;
% 
% function out=get_new(u,para)
% [nn,mm]=size(u);
% m=para.box.m;
% n=para.box.n;
% uu=zeros(nn,mm);
% cal=1;
% for ii=1:m+1
%     for jj=1:n+1
%     id=ii+(jj-1)*(m+1);
%     if jj~=1||jj~=n+1||ii~=1||ii~=m+1
%         uu(cal)=u(id);
%         cal=cal+1;
%     end
%     end
% end
% sum=0;
% for ii=1:nn
%     sum=sum+uu(ii);
% end
% aver=sum/(cal-1);
% out=uu-aver;