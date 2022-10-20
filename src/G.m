function out = G(x,y,para,phi)
g = para.u;
g_ = para.u_;

if phi_func(x,y,para,phi) < 0
	out = g_(x,y);
end

if phi_func(x,y,para,phi) >= 0
	out = g(x,y);
end
