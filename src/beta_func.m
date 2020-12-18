%% Evaluates the coefficient beta at each grid point
% $$ \beta = K \lambda $$
% x and y are the coordinates of the grid point
% The corresponding index locations in the matrix for beta
% are determined in mm and nn respectively.
function out = beta_func(x,y,para,beta)

dx = para.box.dx;
dy = para.box.dy;

%m=para.box.m;
%n=para.box.n;

left = para.box.left;
bottom = para.box.bottom;

mm=round((x-left)/dx)+1;
nn=round((y-bottom)/dy)+1;
 
out=beta(nn,mm);
