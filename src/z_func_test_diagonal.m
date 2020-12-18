%% Specifying the initial position of the water front
% A function describing the initial position of 
% the water front in the shape of a vertical front of a HS cell 
% $$ z(x,y) = y-0.015 $$ 
% This can take array input
%% 
function out = z_func_test_diagonal(x,y)

% square initial profile
 out = max(abs(x-.5),abs(y-.5))-0.2/sqrt(2);