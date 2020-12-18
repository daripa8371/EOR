%% Specifying the initial position of the water front
% A function describing the initial position of 
% the water front in the shape of a circular arc 
% $$ z(x,y) = x^2+y^2-0.015 $$ 
% This can take array input
%% 
function out = z_func_test(x,y)

% out=(x-0.15)^2+(y-0.15)^2 -0.5*(x+y-0.35)^2;
 out=(x).^2+(y).^2-0.015;