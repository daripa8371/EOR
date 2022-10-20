function out = eval_dtconc(dt,x,y,g,u,v,Ds,Dg,Qx,Qy,Gx,Gy)
% function to determine correct dt so that backward differencing along
% characteristics does not produce points outside the domain.
xjump = x-(g*u + Ds*Qx + Dg*Gx)*dt;
yjump = y-(g*v + Ds*Qy + Dg*Gy)*dt;
if xjump < 0 && yjump >= 0 && yjump <= 1
    dtmod = x/(g*u + Ds*Qx + Dg*Gx);
elseif xjump >= 0 && xjump <= 1 && yjump < 0
    dtmod = y/(g*v + Ds*Qy + Dg*Gy);
elseif xjump > 1 && yjump >= 0 && yjump <= 1
    dtmod = (x-1)/(g*u + Ds*Qx + Dg*Gx);
elseif xjump >= 0 && xjump <= 1 && yjump > 1
    dtmod = (y-1)/(g*v + Ds*Qy + Dg*Gy);
elseif xjump < 0 && yjump > 1
    dtmod = min(x/(g*u + Ds*Qx + Dg*Gx),(y-1)/(g*v + Ds*Qy + Dg*Gy));
elseif xjump < 0 && yjump < 0
    dtmod = min(x/(g*u + Ds*Qx + Dg*Gx),y/(g*v + Ds*Qy + Dg*Gy));
elseif xjump > 1 && yjump < 0 
    dtmod = min((x-1)/(g*u + Ds*Qx + Dg*Gx),y/(g*v + Ds*Qy + Dg*Gy));
elseif xjump > 1 && yjump > 1
    dtmod = min((x-1)/(g*u + Ds*Qx + Dg*Gx),(y-1)/(g*v + Ds*Qy + Dg*Gy));
else
    dtmod = dt;
end
    out = dtmod;