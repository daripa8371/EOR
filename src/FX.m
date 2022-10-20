function out=FX(u1,u2,u3,u4,v,miuw,miuo,alpha)
du1=u2-u1;
du2=u3-u2;
du3=u4-u3;

ul=u2+0.5*minmod(du1,du2);
ur=u3-0.5*minmod(du2,du3);
ug=(ul+ur)/2;

ss=ug(:,:,1);
cc=ug(:,:,2);

q1=ss;
q2=ss.*cc;
eigen1= (miuo.*q1.^2)./(miuw.*q1 + miuo.*q1.^3 - 2*miuw.*q1.^2 + miuw.*q1.^3 + alpha.*miuw.*q2 - 2*alpha.*miuw.*q1.*q2 + alpha.*miuw.*q1.^2.*q2);
eigen2=(2*miuo.*miuw.*q1.^3 - 2*miuo.*miuw.*q1.^4 + 2*alpha.*miuo.*miuw.*q1.^2.*q2 - 2*alpha.*miuo.*miuw.*q1.^3.*q2)./(alpha.^2.*miuw.^2.*q1.^4.*q2.^2 - 4*alpha.^2.*miuw.^2.*q1.^3.*q2.^2 + 6*alpha.^2.*miuw.^2.*q1.^2.*q2.^2 - 4*alpha.^2.*miuw.^2.*q1.*q2.^2 + alpha.^2.*miuw.^2.*q2.^2 + 2*alpha.*miuo.*miuw.*q1.^5.*q2 - 4*alpha.*miuo.*miuw.*q1.^4.*q2 + 2*alpha.*miuo.*miuw.*q1.^3.*q2 + 2*alpha.*miuw.^2.*q1.^5.*q2 - 8*alpha.*miuw.^2.*q1.^4.*q2 + 12*alpha.*miuw.^2.*q1.^3.*q2 - 8*alpha.*miuw.^2.*q1.^2.*q2 + 2*alpha.*miuw.^2.*q1.*q2 + miuo.^2.*q1.^6 + 2*miuo.*miuw.*q1.^6 - 4*miuo.*miuw.*q1.^5 + 2*miuo.*miuw.*q1.^4 + miuw.^2.*q1.^6 - 4*miuw.^2.*q1.^5 + 6*miuw.^2.*q1.^4 - 4*miuw.^2.*q1.^3 + miuw.^2.*q1.^2);

lamda=max(abs(v.*eigen1),abs(v.*eigen2));

gl=zeros(size(ul,1),size(ul,2),2);
gr=zeros(size(ul,1),size(ul,2),2);

s=ul(:,:,1);
c=ul(:,:,2);
%gama=ul(:,:,3);

miua=miuw*(1+alpha*c);  % are miua,miuw,miuo matrices? suppose they are all matrices

gl(:,:,1)=v.*s.^2./miua./(s.^2./miua+(1-s).^2/miuo);
gl(:,:,2)=v.*c.*s.^2./miua./(s.^2./miua+(1-s).^2/miuo);

s=ur(:,:,1);
c=ur(:,:,2);
%gama=ur(:,:,3);

miua=miuw*(1+alpha*c);  % are miua,miuw,miuo matrices? suppose they are all matrices

gr(:,:,1)=v.*s.^2./miua./(s.^2./miua+(1-s).^2/miuo);
gr(:,:,2)=v.*c.*s.^2./miua./(s.^2./miua+(1-s).^2/miuo);

out=(gl+gr)/2-cat(3,lamda,lamda)/2.*(ur-ul);
