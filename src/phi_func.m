function out = phi_func(x,y,para,phi)

dx = para.box.dx;
dy = para.box.dy;

left = para.box.left;
bottom = para.box.bottom;

mm=round((x-left)/dx)+1;
nn=round((y-bottom)/dy)+1;
if abs((x-left)/dx+1-mm)>10^(-10)||abs((y-bottom)/dy +1-nn)>10^(-10)

   abs((x-left)/dx+1-mm) 
   abs((y-bottom)/dy +1-nn)
   kkk
end
out=phi(nn,mm);
