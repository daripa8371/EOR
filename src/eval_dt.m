function out = eval_dt(dt,x,y,b_s,u,v)
% function to determine correct dt so that backward differencing along
% characteristics does not produce points outside the domain.

xjump = x-b_s*u*dt;
yjump = y-b_s*v*dt;
if xjump < 0 && yjump >= 0 && yjump <= 1
    dtmod = x/(b_s*u);
elseif xjump >= 0 && xjump <= 1 && yjump < 0
    dtmod = y/(b_s*v);
elseif xjump > 1 && yjump >= 0 && yjump <= 1
    dtmod = (x-1)/(b_s*u);
elseif xjump >= 0 && xjump <= 1 && yjump > 1
    dtmod = (y-1)/(b_s*v);
elseif xjump < 0 && yjump > 1
    dtmod = min(x/(b_s*u),(y-1)/(b_s*v));
elseif xjump < 0 && yjump < 0
    dtmod = min(x/(b_s*u),y/(b_s*v));
elseif xjump > 1 && yjump < 0 
    dtmod = min((x-1)/(b_s*u),y/(b_s*v));
elseif xjump > 1 && yjump > 1
    dtmod = min((x-1)/(b_s*u),(y-1)/(b_s*v));
else
    dtmod = dt;
end
    out = dtmod;
    
