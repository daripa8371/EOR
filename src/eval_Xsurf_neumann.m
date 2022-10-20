function [out1,out2] = eval_Xsurf_neumann(x,y,s,snew,g,f,f_s,D,pc_s,pc_g,u,v,dt,para,flag)
%function to compute shifted coordinates over one time step in the
%characteristic direction
% flag = 1 (b_s)
% flag = 2 (b_c)
% flag = 3 (b_g)

% to check if any dt is missing out in computations 
%dtmod = dt;%999*ones(size(y,1),size(x,2));


if flag ==1
    
    xjump = x-f_s.*u.*dt;
    yjump = y-f_s.*v.*dt;
    
elseif flag == 2
    [sx,sy] = get_gra(s,para);
    [gx,gy] = get_gra(g,para);
    
    xjump = x-((f./snew).*u + (D.*pc_s./snew).*sx + (D.*pc_g./snew).*gx).*dt;
    yjump = y-((f./snew).*v + (D.*pc_s./snew).*sy + (D.*pc_g./snew).*gy).*dt;
    
elseif flag ==3
    [sx,sy] = get_gra(s,para);
    
    xjump = x-((f./snew).*u + (D.*pc_s./snew).*sx).*dt;
    yjump = y-((f./snew).*v + (D.*pc_s./snew).*sy).*dt;
   
end

xmod = x; ymod = y;

for j = 1:size(y,1)
    for i = 1:size(x,2)
        
        if xjump(j,i) <= 1 && yjump(j,i) <= 1
            xmod(j,i) = abs(xjump(j,i));
            ymod(j,i) = abs(yjump(j,i));
        elseif xjump(j,i) > 1 && yjump(j,i) <= 1
            xmod(j,i) = 2 - xjump(j,i);
            ymod(j,i) = abs(yjump(j,i));
        elseif xjump(j,i) <= 1 && yjump(j,i) > 1
            xmod(j,i) = abs(xjump(j,i));
            ymod(j,i) = 2 - yjump(j,i);
        elseif xjump(j,i) > 1 && yjump(j,i) > 1
            xmod(j,i) = 2 - xjump(j,i);
            ymod(j,i) = 2 - yjump(j,i);
        end
        
    end 
end
    
    out1 = xmod;
    out2 = ymod;
    


% if flag ==1
% 
%     for j = 1:size(y,1)
%         for i = 1:size(x,2)
%             dtmod(j,i) = eval_dt(dt(j,i),x(j,i),y(j,i),f_s(j,i),u(j,i),v(j,i));
%         end
%     end
% 
%     xmod = x - f_s.*u.*dtmod;
%     ymod = y - f_s.*v.*dtmod;
%     out1=xmod;
%     out2=ymod;
%     out3=dtmod;
%     
% elseif flag == 2
%     [sx,sy] = get_gra(s,para);
%     [gx,gy] = get_gra(g,para);
%     
%     for j = 1:size(y,1)
%         for i = 1:size(x,2)
%             dtmod(j,i) = eval_dtconc(dt(j,i),x(j,i),y(j,i),f(j,i)/snew(j,i),...
%             u(j,i),v(j,i),D(j,i)*pc_s(j,i)/snew(j,i),D(j,i)*pc_g(j,i)/snew(j,i),sx(j,i),sy(j,i),gx(j,i),gy(j,i));
%         end
%     end
% 
%     
%     xmod = x - (f.*u + D.*pc_s.*sx + D.*pc_g.*gx).*dtmod./snew;
%     ymod = y - (f.*v + D.*pc_s.*sy + D.*pc_g.*gy).*dtmod./snew;
%     out1=xmod;
%     out2=ymod;
%     out3=dtmod;
%     
% elseif flag ==3
%     [sx,sy] = get_gra(s,para);
%     for j = 1:size(y,1)
%         for i = 1:size(x,2)
%             dtmod(j,i) = eval_dtsurf(dt(j,i),x(j,i),y(j,i),f(j,i)/snew(j,i),...
%                 u(j,i),v(j,i),D(j,i)*pc_s(j,i)/snew(j,i),sx(j,i),sy(j,i));
%         end
%     end
% 
%     xmod = x - (f.*u + D.*pc_s.*sx).*dtmod./snew;
%     ymod = y - (f.*v + D.*pc_s.*sy).*dtmod./snew;
%     out1=xmod;
%     out2=ymod;
%     out3=dtmod;
%end