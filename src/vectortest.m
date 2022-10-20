clear all
clc
% a=zeros(5);
% for i = 1:5
%     for j = 1:5
%         a(i,j) = (-1)^(i+j);
%     end
% end
% a
% % b=zeros(10)
% % m=1:10;
% % n=[1:10]';
% % b(m,n) = (-1).^(m+n)
% [r,c] = meshgrid(1:5,1:5);
% b = (-1).^(r+c)
% c = bsxfun(@(A,B) (-1).^(A+B),(1:5),(1:5)')
% isequal(a,c)

I = rand(5); S1=zeros(5);
for i=1:5
    for j=1:5
        if(I(i,j)>=.5||I(i,j)<=0.1)
            S1(i,j) = 1;
        else
            S1(i,j) = 1-I(i,j);
        end
    end
end
S2 = zeros(5); 
%[r,c] = meshgrid(1:5,1:5);
D = (I>=.5)+(I<=0.1);
S2 = D + (~D).*(1-I);
disp(I)
disp(D)
S1
S2
isequal(S1,S2)

% flag = 0;
% 
% D = (phi>10^(-10)||abs(phi)<10^(-10));
%     if flag ==0
%         g0=[];
%         s0=D.*(1-s_0)+(~D);
%         c0=(~D).*c_0;
%     elseif flag == 1
%         %g0=s0;
%         s0=D.*(1-s_0)+(~D);
%         c0=(~D).*c_0;
%         g0=(~D).*g_0;
%     end
% %     if flag == 0
% %         g0 = [];
% %         for i=1:para.box.m+1
% %             for j=1:para.box.n+1
% %                 if phi(j,i)>10^(-10)||abs(phi(j,i))<10^(-10)
% %                     s0(j,i)=1-s_0;
% %                     c0(j,i)=0;
% %                 else
% %                     s0(j,i)=1;
% %                     c0(j,i)=c_0;
% %                 end
% %             end
% %         end
% %     elseif flag == 1    
% %         g0=zeros(para.box.n+1,para.box.m+1);
% %         for i=1:para.box.m+1
% %             for j=1:para.box.n+1
% %                 if phi(j,i)>10^(-10)||abs(phi(j,i))<10^(-10)
% %                     s0(j,i)=1-s_0;
% %                     c0(j,i)=0;
% %                     g0(j,i)=0;
% %                 else
% %                     s0(j,i)=1;
% %                     c0(j,i)=c_0;
% %                     g0(j,i)=g_0;
% %                 end
% %             end
% %         end
% %     end