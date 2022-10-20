clear all
syms q1 q2 miuw alpha miuo u v;
f1=u*q1^2/(miuw*(1+alpha*q2/q1))/(q1^2/(miuw*(1+alpha*q2/q1))+(1-q1^2)/miuo);
f2=u*q2/q1*q1^2/(miuw*(1+alpha*q2/q1))/(q1^2/(miuw*(1+alpha*q2/q1))+(1-q1^2)/miuo);
dfdq=[diff(f1,q1),diff(f1,q2);diff(f2,q1),diff(f2,q2)]
eu=eig(dfdq)