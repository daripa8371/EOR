%% Evaluates the source term f at each grid point
% f is non-zero only at injection and production wells
% x and y are the coordinates of the grid point
% The corresponding index locations in the matrix for f
% are determined in mm and nn respectively.
function out = f_func(x,y,para,f)

dx = para.box.dx;
dy = para.box.dy;

%n=para.box.n;

left = para.box.left;
bottom = para.box.bottom;

mm=round((x-left)/dx)+1;
nn=round((y-bottom)/dy)+1;
 
out=f(nn,mm);