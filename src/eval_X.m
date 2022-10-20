function [out1,out2,out3] = eval_X(x,y,b_s,u,v,dt)
% to check if any dt is missing out in computations 
dtmod = 999*ones(size(y,1),size(x,2));

for j = 1:size(y,1)
    for i = 1:size(x,2)
        dtmod(j,i) = eval_dt(dt(j,i),x(j,i),y(j,i),b_s(j,i),u(j,i),v(j,i));
    end
end
if max(max(dtmod)) == 999
    disp('Error:dtmod = 999')
end

xmod = x - b_s.*u.*dtmod;
ymod = y - b_s.*v.*dtmod;
out1=xmod;
out2=ymod;
out3=dtmod;