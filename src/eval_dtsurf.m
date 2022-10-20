function out = eval_dtsurf(dt,x,y,g,u,v,Ds,Qx,Qy)
% function to determine correct dt so that backward differencing along
% characteristics does not produce points outside the domain.
xjump = x-(g*u + Ds*Qx)*dt;
yjump = y-(g*v + Ds*Qy)*dt;
if xjump < 0 && yjump >= 0 && yjump <= 1
    dtmod = x/(g*u + Ds*Qx);
elseif xjump >= 0 && xjump <= 1 && yjump < 0
    dtmod = y/(g*v + Ds*Qy);
elseif xjump > 1 && yjump >= 0 && yjump <= 1
    dtmod = (x-1)/(g*u + Ds*Qx);
elseif xjump >= 0 && xjump <= 1 && yjump > 1
    dtmod = (y-1)/(g*v + Ds*Qy);
elseif xjump < 0 && yjump > 1
    dtmod = min(x/(g*u + Ds*Qx),(y-1)/(g*v + Ds*Qy));
elseif xjump < 0 && yjump < 0
    dtmod = min(x/(g*u + Ds*Qx),y/(g*v + Ds*Qy));
elseif xjump > 1 && yjump < 0 
    dtmod = min((x-1)/(g*u + Ds*Qx),y/(g*v + Ds*Qy));
elseif xjump > 1 && yjump > 1
    dtmod = min((x-1)/(g*u + Ds*Qx),(y-1)/(g*v + Ds*Qy));
else
    dtmod = dt;
end
    out = dtmod;