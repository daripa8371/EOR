%% Specifying the initial position of the water front
% A function describing the initial position of 
% the water front in the shape of a circular arc 
% $$ z(x,y) = x^2+y^2-0.015 $$ 
% This can take array input
%% 
function out = z_func_test(x,y)


init_front_hs = 0.1;

% out=(x-0.15)^2+(y-0.15)^2 -0.5*(x+y-0.35)^2;



% % perturbed initial saturation front for special fingering simulations
% out = x.^2 + y.^2 - 0.015*(1+0.1*sin(18*atan(y./x)))^2;

% Homogeneous
  out = y - init_front_hs + 0.01*(cos((80*pi*x)));



% Normal unperturbed initial saturation front 
% Quarter 5 spot
%out=(x).^2+(y).^2-0.015;


% % Rectilinear heterogeneous
%     out = y - init_front_hs;