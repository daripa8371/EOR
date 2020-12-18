%%
% Input:
% T is a structure array with fields x & y where
%   T.x contains x coordinates of vertices of an element triangle
%   T.y contains y coordinates of vertices of an element triangle
% beta is the average of the value at the vertices of the 
%   coefficient $$\beta = K(x) \lambda(s,c,\Gamma)$$
% v = ?
% Output: ?
function out = weak1(T,beta,v)

% Computing area of the triangle formed by T
s = polyarea(T.x,T.y); 
% M = matrix where 1st column is x coordinate of element triangle, 2nd column is y coord of elem. triangle and 3rd column is filled with 1's
M = [T.x', T.y', [1;1;1]]; 
M = M^(-1); % finding inverse of the original matrix M (3x3)
M = M([1,2],:); % extracting the first two rows of M to make it 2x3 matrix
vdiff = M*v';  % multiplying M (2x3) with transpose of v (3x1) to get vdiff (2x1)
inte= vdiff'*beta*M*s;
out = [inte,0];
end