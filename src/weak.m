% Input:
% T is a structure array with fields x & y where
%   T.x contains x coordinates of vertices of an element triangle
%   T.y contains y coordinates of vertices of an element triangle
% beta is the average of the value at the vertices of the 
%   coefficient $$\beta = K(x) \lambda(s,c,\Gamma)$$
% v = ?
% Output: ?
function out = weak(T,para,beta,v)
    % evaluating beta at the vertices of the element triangle
    beta1=beta_func(T.x(1),T.y(1),para,beta);

    beta2=beta_func(T.x(2),T.y(2),para,beta);

    beta3=beta_func(T.x(3),T.y(3),para,beta);
    
    % computing average of the beta values at the vertices
    beta10=(beta1+beta2+beta3)/3;
    out = weak1(T,beta10,v);
   
end