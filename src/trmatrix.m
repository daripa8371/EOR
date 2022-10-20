function out = trmatrix(T,f1,f2,f3,v)
s = polyarea(T.x,T.y);
fc=(f1+f2+f3)/3;
vc = (v(1)+v(2)+v(3))/3;
% out = (f1*v(1)+f2*v(2)+f3*v(3))*s/3;
f4 = (f2+f3)/2;
f5 = (f1+f3)/2;
f6 = (f1+f2)/2;



v4 = (v(2)+v(3))/2;
v5 = (v(3)+v(1))/2;
v6 = (v(1)+v(2))/2;

out = (f4*v4+f5*v5+f6*v6+fc*vc)*s/4;
