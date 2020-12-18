%% Specifying the initial position of the water front
% A function describing the initial position of 
% the water front in the shape of a circle at the center of a square  
% $$ z(x,y) = (x-.5)^2+ (y-0.5)^2 - 0.015 $$ 
% This can take array input
%% 
function out = z_func_test_parallel(x,y)
% circular initial profile
 %out=(x-0.5).^2+(y-0.5).^2-0.025;
 % square initial profile
 out = abs(x-.5) + abs(y-.5) - 0.2;
 