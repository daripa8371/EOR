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
M = [T.x', T.y', [1;1;1]];
M = M^(-1);
M = M([1,2],:);
vdiff = M*v';
inte= vdiff'*beta*M*s;
out = [inte,0];
end