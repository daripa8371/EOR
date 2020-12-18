%% Specifying the initial position of the water front
% A function describing the initial position of 
% the water front in the shape of a vertical front of a HS cell 
% $$ z(x,y) = x - 0.1 $$ for unperturbed interface
% $$ z(x,y) = x - 0.1 + 0.01*(cos((16*pi*y))) $$ for perturbed interface
% This can take array input
%% 
function out = z_func_test_hs(x,y)
global init_front_hs


 out = x - init_front_hs + 0.01*(cos((16*pi*y)));