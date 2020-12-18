%% Setting up triangulations for the FEM grid
% U = cell array with each element = array of vertices of Upper Triangle of
% the rectangular cell 
% L = cell array with each element = array of vertices of Lower Triangle of
% the rectangular cell
% At every point (i,j), U{i,j} & L{i,j} are cells with coordinates of vertices
% of the two triangles obtained by bisecting the rectangle starting at
% (i,j). The bisection line goes from NW to SE.
function [U,L] = setTri(para)
m = para.box.m;
n = para.box.n;
dx = para.box.dx;
dy = para.box.dy;

left = para.box.left;
bottom = para.box.bottom;


U = cell(m,n);
L = cell(m,n);
for j = 1:m
    for k = 1:n
        x1 = left + (j-1)*dx;
        y1 = bottom + (k-1)*dy;
        x2 = left + j*dx;
        y2 = y1;
        x3 = x1;
        y3 = bottom + k*dy;
        x4 = x2;
        y4 = y3;
               
        l.x = [x1,x2,x3];
        l.y = [y1,y2,y3];
        u.x = [x4,x3,x2];
        u.y = [y4,y3,y2];
        

        %-------------------------------------
        U{j,k} = u;
        L{j,k} = l;
        clear u;
        clear l;
    end
end











