function out=get_pc(KK)
porousity=0.1;
tensor=0.1;
J=0.1;
out=(porousity./KK).^(1/2).*tensor.*cos(pi/4).*J;